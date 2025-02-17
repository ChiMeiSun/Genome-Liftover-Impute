# snakemake --dag  | dot -Tsvg > dag.svg
# snakemake --rulegraph  | dot -Tsvg > rulegraph.svg
configfile: "config.yaml"

rule all:
    input:
        (config["Output"]["target_impute"] + ".vcf.gz"),
        config["Output"]["genolift"] + ".vcf.gz",
        expand(config["Impacc"]["dir"]+"/CV_rep{rep}_fold{fold}.pdf", rep=config["Impacc"]["rep"], fold=config["Impacc"]["fold"])

rule make_chain:
    input:
        old = config["Assembly"]["old"],
        new = config["Assembly"]["new"],
        parallel = config["Packages"]["parallel"],
    output:
        chain = protected(config["Output"]["chain"] + ".chain")
    log: "Logs/make_chain.log"  
    benchmark: "Logs/make_chain.bench" 
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["chainmem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        fasplitsize = config["Params"]["fasplitsize"]
    shadow: "shallow"
    shell:
        """
        ( tar -xjvf {input.parallel} -C .
            faToTwoBit {input.old} OldAsby.2bit
            faToTwoBit {input.new} NewAsby.2bit

            twoBitInfo OldAsby.2bit OldAsby.chrom.sizes
            twoBitInfo NewAsby.2bit NewAsby.chrom.sizes

            mkdir fasplit
            faSplit size {input.new} {params.fasplitsize} fasplit/size -lift=fasplit/lift.lft

            mkdir psl
            ls fasplit/*.fa | parallel --jobs 0 --use-cpus-instead-of-cores --progress \
            'pblat OldAsby.2bit {{}} psl/{{/.}}.psl \
            -tileSize=12 -minIdentity=98 -noHead -minScore=100 -threads=20'
            # optional: --use-cpus-instead-of-cores, --progress(print progress)

            for i in psl/*.psl; \
            do liftUp -pslQ psl/`basename $i .psl`.lft.psl fasplit/lift.lft warn $i; done

            mkdir chain
            for i in psl/*.lft.psl; \
            do axtChain -linearGap=medium \
            -psl $i OldAsby.2bit NewAsby.2bit chain/`basename $i .lft.psl`.chain; done

            mkdir chainMerge
            chainMergeSort chain/*.chain | chainSplit chainMerge stdin -lump=1000
            cat chainMerge/*.chain > all.chain
            chainSort all.chain all.sorted.chain
            # stdin: standard input

            mkdir net
            chainNet all.sorted.chain OldAsby.chrom.sizes NewAsby.chrom.sizes net/all.net /dev/null
            # /dev/null: discard query.net (the output over the query genome)

            netChainSubset net/all.net all.chain {output.chain}
        ) &> {log}
        """


rule bcftools_liftover:
    input:
        bvcf = config["Data"]["liftgt"],
        chain = (config["Output"]["chain"] + ".chain"),
        oldf = config["Assembly"]["old"],
        newf = config["Assembly"]["new"],
        bcftools = config["Packages"]["bcftools"],
    output:
        lift = protected(config["Output"]["genolift"] + ".vcf.gz"),
        rej = protected(config["Output"]["genorej"] + ".vcf.gz"),
        csi = protected(config["Output"]["genolift"] + ".vcf.gz.csi"),
    log: "Logs/bcftools_liftover.log"  
    benchmark: "Logs/bcftools_liftover.bench" 
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["mem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    shell:
        """
        ( 
          if [ ! -d bcftools-1.20 ]; then
            tar -xjvf {input.bcftools} -C .
          fi
        
        #bcftools-1.20/bcftools +bcftools-1.20/plugin_score/score.so

            grep -v ^## {input.chain} > for_bcftools.chain
            echo for_bcftools.chain done
            
            FILETYPE=""
            if [[ "{input.bvcf}" == *.bcf ]]; then
                FILETYPE="-Ou"
            elif [[ "{input.bvcf}" == *.bcf.gz ]]; then
                FILETYPE="-Ob"
            elif [[ "{input.bvcf}" == *.vcf.gz ]]; then
                FILETYPE="-Oz"
            elif [[ "{input.bvcf}" == *.vcf ]]; then
                FILETYPE="-Ov"
            fi

            bcftools-1.20/bcftools +bcftools-1.20/plugin_score/liftover.so $FILETYPE {input.bvcf} \
            -- -c for_bcftools.chain \
            -s {input.oldf} \
            -f {input.newf} \
            --lift-mt --write-src --write-fail \
            --reject {output.rej} --reject-type z \
            | bcftools-1.20/bcftools sort \
            -Oz -o {output.lift} -W

            rm for_bcftools.chain
        ) &> {log}
        """
# https://software.broadinstitute.org/software/score/


rule beagle_phase_ref:
    input: 
        vcf = config["Data"]["ref"],
        jar = config["Packages"]["beagle"],
    output:
        contig = protected(config["Output"]["ref_phase"] + "_contiglist.txt"),
        phased = protected(config["Output"]["ref_phase"] + ".vcf.gz"),
        csi = protected(config["Output"]["ref_phase"] + ".vcf.gz.csi"),
    log: "Logs/beagle_phase_ref.log"
    benchmark: "Logs/beagle_phase_ref.bench"
    shadow:"shallow"
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["mem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        mincount = config["Params"]["mincount"],
        phased = config["Output"]["ref_phase"]
    shell:
        """
        set +euo pipefail
        ( ## remove contigs with too little variants, remove multialleleic snps 
            if [[ -f "{input.vcf}.csi" || -f "{input.vcf}.tbi" ]]; then
                rm -f {input.vcf}.csi {input.vcf}.tbi
            fi
            bcftools-1.20/bcftools index -f {input.vcf}

            bcftools-1.20/bcftools query -f '%CHROM\n' {input.vcf} | uniq -c | awk '$1>{params.mincount} {{print $0}}' \
            > {output.contig}
            
            keepcontig=$(awk '{{print $2}}' {output.contig} | paste -sd ',')
            echo $keepcontig   

            bcftools-1.20/bcftools view -r $keepcontig {input.vcf} | \
            bcftools-1.20/bcftools view -M2 -m2 | \
            bcftools-1.20/bcftools norm -d all | \
            bcftools-1.20/bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' \
            -Oz -o tmp_qc.vcf.gz --threads {threads}

            java -Xmx$(( {resources.mem_mb} / 1000 - 2 ))g -XX:ConcGCThreads={threads} -jar {input.jar} \
            gt=tmp_qc.vcf.gz out={params.phased}

            bcftools-1.20/bcftools index -f {output.phased}
            
        ) &> {log}
        """


rule beagle_phase_gt:
    input: 
        vcf = config["Data"]["target"],
        jar = config["Packages"]["beagle"],
    output:
        contig = protected(config["Output"]["target_phase"] + "_contiglist.txt"),
        phased = protected(config["Output"]["target_phase"] + ".vcf.gz"),
        csi = protected(config["Output"]["target_phase"] + ".vcf.gz.csi"),
    log: "Logs/beagle_phase_gt.log"
    benchmark: "Logs/beagle_phase_gt.bench"
    shadow:"shallow"
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["mem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        mincount = config["Params"]["mincount"],
        phased = config["Output"]["target_phase"]
    shell:
        """
        set +euo pipefail
        ( ## remove contigs with too little variants, remove multialleleic snps 
            if [[ -f "{input.vcf}.csi" || -f "{input.vcf}.tbi" ]]; then
                rm -f {input.vcf}.csi {input.vcf}.tbi
            fi
            bcftools-1.20/bcftools index -f {input.vcf}

            bcftools-1.20/bcftools query -f '%CHROM\n' {input.vcf} | uniq -c | awk '$1>{params.mincount} {{print $0}}' \
            > {output.contig}
            
            keepcontig=$(awk '{{print $2}}' {output.contig} | paste -sd ',')
            echo $keepcontig   

            bcftools-1.20/bcftools view -r $keepcontig {input.vcf} | \
            bcftools-1.20/bcftools view -M2 -m2 | \
            bcftools-1.20/bcftools norm -d all | \
            bcftools-1.20/bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' \
            -Oz -o tmp_qc.vcf.gz --threads {threads}

            java -Xmx$(( {resources.mem_mb} / 1000 - 2 ))g -XX:ConcGCThreads={threads} -jar {input.jar} \
            gt=tmp_qc.vcf.gz out={params.phased}

            bcftools-1.20/bcftools index -f {output.phased}
        ) &> {log} || true
        """

rule check_preimpute:
    input:
        contig_ref = (config["Output"]["ref_phase"] + "_contiglist.txt"),
        contig_gt = (config["Output"]["target_phase"] + "_contiglist.txt"),
        jar = config["Packages"]["cfgt"],
        ref = config["Output"]["ref_phase"] + ".vcf.gz",
        gt = config["Output"]["target_phase"] + ".vcf.gz",
    output:
        out = config["Output"]["target_phase"] + "_aligned.vcf.gz",
    log: "Logs/check_preimpute.log"
    benchmark: "Logs/check_preimpute.bench"
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["mem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        path = "Output/Phase/conformgt/",
    shell:
        """
        set +euo pipefail
        ( ## ref must include all contigs as in gt, otherwise remove those not from gt
            comchr=$(comm -12 <(awk '{{print$2}}' {input.contig_ref} | sort) <(awk '{{print $2}}' {input.contig_gt} | sort) | sort -n | paste -sd " ")
            echo $comchr

        if [ -d {params.path} ]; then
            rm -r {params.path} 
        fi
        mkdir {params.path}

        for i in $comchr; do
            echo $i
            java -jar {input.jar} chrom=$i match=POS \
                ref={input.ref} \
                gt={input.gt} \
                out={params.path}$i

            bcftools-1.20/bcftools index -f {params.path}$i.vcf.gz
        done
        
        files=$(ls {params.path}*.vcf.gz | sort -n)
        echo $files
        bcftools-1.20/bcftools concat $files -Oz -o {output.out}

        ) &> {log}
        """
# https://faculty.washington.edu/browning/conform-gt.html


rule beagle_impute:
    input: 
        bref = config["Packages"]["bref"],
        ref = (config["Output"]["ref_phase"] + ".vcf.gz"),
        jar = config["Packages"]["beagle"],
        gt = (config["Output"]["target_phase"] + "_aligned.vcf.gz")
    output:
        bref3 = protected(config["Output"]["ref_phase"] + ".bref3"),
        imputed = protected(config["Output"]["target_impute"] + ".vcf.gz"),
        csi = protected(config["Output"]["target_impute"] + ".vcf.gz.csi")
    log: "Logs/beagle_impute.log"
    benchmark: "Logs/beagle_impute.bench"
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["mem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        imputed = config["Output"]["target_impute"],
        tf = config["Params"]["impute"]
    shell:
        """
        set +euo pipefail
        (   java -jar {input.bref} {input.ref} > {output.bref3}
            echo convert vcf to bref3 done

            java -Xmx$(( {resources.mem_mb} / 1000 - 2 ))g -XX:ConcGCThreads={threads} -jar {input.jar} \
            impute={params.tf} \
            gt={input.gt} \
            ref={output.bref3} \
            out={params.imputed}

            bcftools-1.20/bcftools index -f {output.imputed}
            ) &> {log}
        """

rule ImputationAccuracy:
    input:
        script = config["Packages"]["cv"],
        target = (config["Output"]["target_phase"] + "_aligned.vcf.gz"),
        ref = (config["Output"]["ref_phase"] + ".bref3"),
        jar = config["Packages"]["beagle"],
        parallel = config["Packages"]["parallel"],
    output:
        plot = protected(config["Impacc"]["dir"]+"/CV_rep{rep}_fold{fold}.pdf"),
        res = protected(config["Impacc"]["dir"]+"/CV_rep{rep}_fold{fold}.txt"),
    log: ("Logs/"+config["Impacc"]["dir"]+"/CV_rep{rep}_fold{fold}.log")
    benchmark: ("Logs/"+config["Impacc"]["dir"]+"/CV_rep{rep}_fold{fold}.bench")
    threads: config["Resources"]["threads"]
    resources:
        time = lambda wildcards, attempt: int(attempt) *config["Resources"]["time"],
        mem_mb = lambda wildcards, attempt: int(attempt) *config["Resources"]["mem"],
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    shell:
        """
        set +euo pipefail
        ( Rscript --vanilla {input.script} {input.target} {input.ref} {input.jar} \
        {wildcards.rep} {wildcards.fold} {resources.mem_mb} {threads} {output.plot} {output.res} {input.parallel}
        ) &> {log} 
        """ 
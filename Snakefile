

rule all:
    input:
        # "Output/Imputed/bcftools_lifted_imputed_gt.vcf.gz"

rule make_chain:
    input:
        old = "Data/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.primary_assembly.25.fa",
        new = "Data/Gallus_gallus.GRCg6a.dna_rm.chromosome.25.fa",
        parallel = "Packages/parallel-latest.tar.bz2"
    output:
        chain = protected("Output/Chain/GRCg7bChr25_to_GRCg6aChr25.chain")
    log: "Logs/make_chain.log"  
    benchmark: "Logs/make_chain.bench" 
    conda: "Chain_Liftover"
    threads: 20
    resources:
        time = lambda wildcards, attempt: int(attempt) *24*2* 60,
        mem_mb = lambda wildcards, attempt: int(attempt) *1000*1000,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        oldname = "GRCg7bChr25",
        newname = "GRCg6aChr25",
        fasplitsize = "10000000"
    shadow: "shallow"
    shell:
        """
        ( tar -xjvf {input.parallel} -C .
            faToTwoBit {input.old} {params.oldname}.2bit
            faToTwoBit {input.new} {params.newname}.2bit

            twoBitInfo {params.oldname}.2bit {params.oldname}.chrom.sizes
            twoBitInfo {params.newname}.2bit {params.newname}.chrom.sizes

            mkdir fasplit
            faSplit size {input.new} {params.fasplitsize} fasplit/size -lift=fasplit/lift.lft

            mkdir psl
            ls fasplit/*.fa | parallel --jobs 0 --use-cpus-instead-of-cores --progress \
            'pblat {params.oldname}.2bit {{}} psl/{{/.}}.psl \
            -tileSize=12 -minIdentity=98 -noHead -minScore=100 -threads=20'
            # optional: --use-cpus-instead-of-cores, --progress(print progress)

            for i in psl/*.psl; \
            do liftUp -pslQ psl/`basename $i .psl`.lft.psl fasplit/lift.lft warn $i; done

            mkdir chain
            for i in psl/*.lft.psl; \
            do axtChain -linearGap=medium \
            -psl $i {params.oldname}.2bit {params.newname}.2bit chain/`basename $i .lft.psl`.chain; done

            mkdir chainMerge
            chainMergeSort chain/*.chain | chainSplit chainMerge stdin -lump=1000
            cat chainMerge/*.chain > all.chain
            chainSort all.chain all.sorted.chain
            # stdin: standard input

            mkdir net
            chainNet all.sorted.chain {params.oldname}.chrom.sizes {params.newname}.chrom.sizes net/all.net /dev/null
            # /dev/null: discard query.net (the output over the query genome)

            netChainSubset net/all.net all.chain {output.chain}
        ) &> {log}
        """


rule bcftools_liftover:
    input:
        bvcf = "Data/autosomes.ARaucanaWL.annot.merged.filtered.bcf",
        chain = "Output/Chain/GRCg7bChr25_to_GRCg6aChr25.chain",
        oldf = "Data/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.primary_assembly.25.fa",
        newf = "Data/Gallus_gallus.GRCg6a.dna_rm.chromosome.25.fa",
        bcftools = "Packages/bcftools-1.20.tar.bz2",
    output:
        bcfchain = "Output/Chain/GRCg7b_to_GRCg6a_for_bcftools.chain",
        lift = protected("Output/Liftover/bcftools_lifted.vcf.gz"),
        rej = protected("Output/Liftover/bcftools_rejected_variants.vcf.gz"),
        csi = protected("Output/Liftover/bcftools_lifted.vcf.gz.csi"),
    log: "Logs/bcftools_liftover.log"  
    benchmark: "Logs/bcftools_liftover.bench" 
    conda: "Chain_Liftover"
    threads: 20
    resources:
        mem_mb = lambda wildcards, attempt: int(attempt) * 30 * 1000,
        time = lambda wildcards, attempt: int(attempt) *6* 60,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    shell:
        """
        ( 
          if [ ! -d bcftools-1.20 ]; then
            tar -xjvf {input.bcftools} -C .
          fi
        
        #bcftools-1.20/bcftools +bcftools-1.20/plugin_score/score.so

            grep -v ^## {input.chain} > {output.bcfchain}
            echo bcftools_lift_chain done
            
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
            -- -c {output.bcfchain} \
            -s {input.oldf} \
            -f {input.newf} \
            --lift-mt --write-src --write-fail \
            --reject {output.rej} --reject-type z \
            | bcftools-1.20/bcftools sort \
            -Oz -o {output.lift} -W
        ) &> {log}
        """
# https://software.broadinstitute.org/software/score/


rule beagle_phase_ref:
    input: 
        vcf = "Output/Liftover/bcftools_lifted.vcf.gz",
        jar = "Packages/beagle.22Jul22.46e.jar"
    output:
        contig = protected("Output/Liftover/bcftools_lifted_contiglist.txt"),
        phased = protected("Output/Phased/bcftools_lifted_qc_phased.vcf.gz"),
        csi = protected("Output/Phased/bcftools_lifted_qc_phased.vcf.gz.csi"),
    log: "Logs/beagle_phase_ref.log"
    benchmark: "Logs/beagle_phase_ref.bench"
    shadow:"shallow"
    conda: "Chain_Liftover"
    threads: 20
    resources:
        mem_mb = lambda wildcards, attempt: int(attempt) *100* 1000,
        time = lambda wildcards, attempt: int(attempt) *3* 60,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        mincount = 1000,
        phased = "Output/Phased/bcftools_lifted_qc_phased"
    shell:
        """
        set +o pipefail
        ( ## remove contigs with too little variants, remove multialleleic snps 
        
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


rule beagle_phase_gt:
    input: 
        vcf = "Data/IMAGE.vcf",
        jar = "Packages/beagle.22Jul22.46e.jar"
    output:
        out = "Output/Phased/IMAGE_phased.vcf.gz",
        csi = "Output/Phased/IMAGE_phased.vcf.gz.csi",
        contig = "Output/Phased/IMAGE_contiglist.txt",
    log: "Logs/beagle_phase_gt.log"
    benchmark: "Logs/beagle_phase_gt.bench"
    conda: "Chain_Liftover"
    threads: 20
    resources:
        mem_mb = lambda wildcards, attempt: int(attempt) * 20*1000,
        time = lambda wildcards, attempt: int(attempt) *3* 60,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        out = "Output/Phased/IMAGE_phased"
    shell:
        """
        (   java -Xmx$(( {resources.mem_mb} / 1000 - 2 ))g -XX:ConcGCThreads={threads} -jar {input.jar} \
            nthreads={threads} gt={input.vcf} out={params.out}

            bcftools-1.20/bcftools index -f {output.out}
            bcftools-1.20/bcftools query -f '%CHROM\n' {output.out} | uniq > {output.contig}
            ) &> {log}
        """

rule check_preimpute:
    input:
        contig_ref = "Output/Liftover/bcftools_lifted_contiglist.txt",
        contig_gt = "Output/Phased/IMAGE_contiglist.txt",
        jar = "Packages/conform-gt.24May16.cee.jar",
        ref = "Output/Phased/bcftools_lifted_qc_phased.vcf.gz",
        gt = "Output/Phased/IMAGE_phased.vcf.gz",
    output:
        out = protected("Output/Phased/conformgt/IMAGE_phased_aligned.vcf.gz"),
        gtc = "Output/Phased/bcfgtcheck_seq_IMAGE_phased_arau.txt"
    log: "Logs/check_conformgt.log"
    benchmark: "Logs/check_conformgt.bench"
    conda: "Chain_Liftover"
    threads: 20
    resources:
        mem_mb = lambda wildcards, attempt: int(attempt) *100* 1000,
        time = lambda wildcards, attempt: int(attempt) * 2*24*60,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        path = "Output/Phased/conformgt/",
        name = "IMAGE_phased_aligned.vcf.gz"
    shell:
        """
        set +o pipefail
        ( ## ref must include all contigs as in gt, otherwise remove those not from gt
            comchr=$(comm -12 <(awk '{{print$2}}' {input.contig_ref} | sort) <(sort {input.contig_gt}) | sort -n | paste -sd " ")
            echo $comchr

            samref=$(bcftools-1.20/bcftools query -l {input.ref} | head -6 | paste -sd ',')        
            echo $samref
            samgt=$(bcftools-1.20/bcftools query -l {input.gt} | head -6 | paste -sd ',')
            echo $samgt

        rm -r {params.path}
        mkdir {params.path}
        for i in $comchr; do
            echo $i
            java -jar {input.jar} chrom=$i match=POS \
                ref={input.ref} \
                gt={input.gt} \
                out={params.path}$i

            bcftools-1.20/bcftools index -f {params.path}$i.vcf.gz
        done

            bcftools-1.20/bcftools gtcheck -E 0 -u GT,GT \
            -s gt:$samgt -s qry:$samref \
            -g {input.gt} {input.ref} > {output.gtc}
        
        cd {params.path}
        files=$(ls *.vcf.gz | sort -n)
        echo $files
        bcftools-1.20/bcftools concat $files -Oz -o {params.name}

        ) &> {log} || true
        """

rule beagle_impute:
    input: 
        bref = "Packages/bref3.27May24.118.jar",
        ref = "Output/Phased/bcftools_lifted_qc_phased.vcf.gz",
        jar = "Packages/beagle.22Jul22.46e.jar",
        gt = "Output/Phased/conformgt/IMAGE_phased_aligned.vcf.gz"
    output:
        bref3 = protected("Output/Phased/bcftools_lifted_qc_phased.bref3"),
        imputed = protected("Output/Imputed/bcftools_lifted_imputed_gt.vcf.gz"),
        csi = protected("Output/Imputed/bcftools_lifted_imputed_gt.vcf.gz.csi")
    log: "Logs/beagle_impute.log"
    benchmark: "Logs/beagle_impute.bench"
    conda: "Chain_Liftover"
    threads: 20
    resources:
        mem_mb = lambda wildcards, attempt: int(attempt) *100* 1000,
        time = lambda wildcards, attempt: int(attempt) *2* 60,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    params:
        imputed="Output/Imputed/bcftools_lifted_imputed_gt"
    shell:
        """
        (
            java -jar {input.bref} {input.ref} > {output.bref3}
            echo convert vcf to bref3 done

            java -Xmx$(( {resources.mem_mb} / 1000 - 2 ))g -XX:ConcGCThreads={threads} -jar {input.jar} \
            impute=true \
            gt={input.gt} \
            ref={output.bref3} \
            out={params.imputed}

            bcftools-1.20/bcftools index -f {output.imputed}
            ) &> {log}
        """

rule check_imputation_CV:
    input:
        script = "Packages/crossVal_imputation.r",
        gt = "Output/Phased/conformgt/IMAGE_phased_aligned.vcf.gz",
        ref = "Output/Phased/bcftools_lifted_qc_phased.bref3",
        jar = "Packages/beagle.22Jul22.46e.jar" 
    output:
        plot = protected("../RESULTS/check_imputation/CV_gt_imp_rep{rep}_fold{fold}.pdf"),
        res = protected("../RESULTS/check_imputation/CV_gt_imp_rep{rep}_fold{fold}.txt")
    log: "../LOGS/check_imputation/check_imputation_CV_rep{rep}_fold{fold}.log"
    benchmark: "../LOGS/check_imputation/check_imputation_CV_rep{rep}_fold{fold}.bench"
    threads: 30
    conda: "IMAGE_SCM"   
    resources:
        mem_mb = lambda wildcards, attempt: int(attempt) *200* 1000,
        time = lambda wildcards, attempt: int(attempt) * 3*24*60,
        temp = lambda wildcards,attempt: "${TMP_LOCAL}" if int(attempt)==1 else 'TEMP/'
    shell:
        """
        ( Rscript --vanilla {input.script} {input.gt} {input.ref} {input.jar} \
        {wildcards.rep} {wildcards.fold} {resources.mem_mb} {threads} {output.plot} {output.res}
        ) &> {log} 
        """ 
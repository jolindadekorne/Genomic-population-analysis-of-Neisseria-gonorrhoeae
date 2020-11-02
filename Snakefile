IDS, = glob_wildcards("raw_data/{id}_R1.fastq.gz")

rule all:
        input:
                expand("fastp_out/{sample}_fastp.json", sample = IDS),
                "multiqc_fastp_out",
                expand("skesa_out/{sample}.fasta", sample = IDS),
                expand("quast_out/{sample}", sample=IDS),
                "multiqc_quast_out",
                expand("sorted_bam_FA1090/{sample}.sorted.bam", sample = IDS),
                expand("sorted_bam_FA1090/samtools_coverage_{sample}.txt", sample = IDS),
                "mlst.tsv",
		"ngmast.tsv",
		expand("ariba/{sample}", sample=IDS),
		expand("snippy_refFA1090_out/{sample}", sample=IDS),
                "snippy_refFA1090_out/clean.full.aln",
		"gubbins_snippyrefFA1090_out/gubbins_snippyrefFA1090.node_labelled.final_tree.tre"

rule fastp:
        input:
                fw = "raw_data/{sample}_R1.fastq.gz",
                rv = "raw_data/{sample}_R2.fastq.gz"
        output:
                fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
                rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz",
                json = "fastp_out/{sample}_fastp.json",
                html = "fastp_out/{sample}_fastp.html"
        params:
                general = "--disable_length_filtering",
                compression_level = 9
        log:
                "logs/fastp/fastp_{sample}.log"
        threads: 4
        shell:
                """
                fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}

                """

rule multiqc_fastp:
        input:
                expand("fastp_out/{sample}_fastp.json", sample=IDS)
        output:
                directory("multiqc_fastp_out")
        log:
                "logs/multiqc/multiqc_fastp.log"
        shell:
                """
                mkdir -p {output}
                multiqc {input} --outdir {output} 2>&1>{log}
                """

rule skesa:
        input:
                fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
                rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz"
        output:
                assembly = "skesa_out/{sample}.fasta"
        params:
                min_length = "500"
        log:
                "logs/skesa/skesa_{sample}.log"
        threads: 6
        shell:
                """
                skesa --fastq {input.fw},{input.rv} --use_paired_ends --min_contig {params.min_length} 1> {output.assembly} 2>{log}
                """

rule quast:
        input:
                assembly = "skesa_out/{sample}.fasta"
        output:
                directory("quast_out/{sample}")
        log:
                "logs/quast/quast_{sample}.log"
        shell:
                """
                quast.py -o {output} {input.assembly}
                """

rule multiqc_quast:
        input:
                expand("quast_out/{sample}/report.tsv", sample=IDS)
        output:
		directory("multiqc_quast_out")
        log:
                "logs/multiqc_quast/multiqc_quast.log"
        shell:
                """
                multiqc {input} -o {output} 2>&1>{log}
                """

rule samtools:
        input:
                fw = "raw_data/{sample}_R1.fastq.gz",
                rv = "raw_data/{sample}_R2.fastq.gz",
                ref = "NgRefFA1090.fna"
        output:
                sort_bam = "sorted_bam_FA1090/{sample}.sorted.bam"
        log:
                "logs/samtools/samtools_{sample}.log"
        params:
                minimap = "sr"
        threads: 6
        shell:
                """
                minimap2 -ax {params.minimap} -t {threads} {input.ref} {input.fw} {input.rv} | samtools view -S -b | samtools sort --threads {threads} -o {output.sort_bam} 1>{output} 2>{log}
                """

rule samtools_coverage:
        input:
                "sorted_bam_FA1090/{sample}.sorted.bam"
        output:
                filelist = "sorted_bam_FA1090/{sample}.txt",
                cov = "sorted_bam_FA1090/samtools_coverage_{sample}.txt"
        log:
                "logs/samtools_coverage_{sample}.log"
        shell:
                """
                ls {input} > {output.filelist}
                samtools coverage -b {output.filelist} -o {output.cov}
                """

rule mlst:
        input:
                expand("skesa_out/{sample}.fasta", sample=IDS)
        output:
                file = "mlst.tsv"
        params:
                scheme = "neisseria"
        log:
                "logs/mlst/mlst.txt"
        shell:
                """
                mlst --legacy --scheme {params.scheme} {input} 1>{output.file} 2>{log}
                """

rule ngmast:
	input:
                expand("skesa_out/{sample}.fasta", sample=IDS)
        output:
                file = "ngmast.tsv"
        log:
                "logs/ngmast/ngmast.txt"
        shell:
                """
                ngmaster {input} 1>{output.file} 2>{log}
                """

rule ariba:
	input:
		fw = "raw_data/{sample}_R1.fastq.gz",
                rv = "raw_data/{sample}_R2.fastq.gz"	
	output:
		directory("ariba/{sample}")
	params:
		refseq="ariba/23S_seq.prepareref"
	shell:
		"""
		ariba run {params.refseq} {input.fw} {input.rv} {output}
		"""

rule snippy_FA1090:
        input:
                fw = "trimmed_illumina/{sample}_fastpout_1.fastq.gz",
                rv = "trimmed_illumina/{sample}_fastpout_2.fastq.gz"
        output:
                directory("snippy_refFA1090_out/{sample}")
        params:
                outdir = directory("snippy_refFA1090_out"),
                ref = "NgRefFA1090.gbk"
        log:    "logs/snippy_refFA1090/snippy_refFA1090_{sample}.log"
        threads: 16
        shell:
                """
                mkdir -p {params.outdir}
                snippy --cpus {threads} --ref {params.ref} --R1 {input.fw} --R2 {input.rv} --outdir {output} 2>{log}
                """

rule snippy_aln:
        input:
                corein = expand("snippy_refFA1090_out/{sample}", sample = IDS)
        output:
                cleanout = "snippy_refFA1090_out/clean.full.aln"
        params:
                outdir = directory("snippy_refFA1090_out"),
                ref = "NgRefFA1090.gbk",
                prefix = "core",
                cleanin = "core.full.aln"
        shell:
                """
                snippy-core {input.corein} --ref {params.ref} --prefix {params.prefix}
		snippy-clean_full_aln {params.cleanin} > {output.cleanout}
                """

rule gubbins:
	input:
		gubbins = "snippy_refFA1090_out/clean.full.aln"
	output:
		"gubbins_snippyrefFA1090_out/gubbins_snippyrefFA1090.node_labelled.final_tree.tre"
	params:	
		model = "GTRGAMMA",
		prefix = "gubbins_snippyrefFA1090"
	log:	"logs/gubbins/gubbins_snippyrefFA1090.log"
	shell:
		"""
		run_gubbins.py {input.gubbins} --prefix {params.prefix} --raxml_model {params.model} 2>{log}
		"""



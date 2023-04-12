#! /usr/bin/env nextflow

// // Convert relative path to absolute path
// if (params.primerBedFile[0]=='/' || params.primerBedFile[0]=='~')
// 	primerBedFile=params.primerBedFile
// else
// 	primerBedFile="$launchDir/$params.primerBedFile"


// params.str = 'Hello world!'

// process splitLetters {
//     output:
//         path 'chunk_*'

//   """
//   printf '${params.str}' | split -b 6 - chunk_
//   """
// }

// process convertToUpper {
//     input:
//         path x
//     output:
//         stdout

//     """
//     cat $x | tr '[a-z]' '[A-Z]'
//     """
// }

// process foo {
//     publishDir '/data/chunks'

//     output:
//     path 'chunk_*'

//     '''
//     printf 'Hola' | split -b 1 - chunk_
//     '''
// }

// process bar {
//     output:
//         stdout

//     l = Channel.of("doing","something","1","2","3","4")

//     """
//     echo l
//     """
// }

process cpu_example {
  cpus 8
  executor 'sge'

  """
  blastp -query input_sequence -num_threads ${task.cpus}
  """
}

// ^^^ tests/examples ^^^ vvv real stuff vvv //

def printUsage () {
	println("Example usage: ")
	println("\tnextflow run ./analyzeReads.nf --plate plate_name --fastqs path/to/fastq_dir --meta path/to/metadata.csv --outdir path/to/outputDir --scheme path/to/primer_scheme_dir")
	println()
	println("Run the artic pipeline and variant analyses.")
	println()
}

process check4kraken {
    cpus 16
    output:
        env kraken_installed
    
    shell:
    '''
    if [[ ! -d !{params.kraken_db} ]]; then
    kraken2-build --db !{params.kraken_db} --download-taxonomy
    kraken2-build --db !{params.kraken_db} --download-library human
    kraken2-build --db !{params.kraken_db} --build --threads !{task.cpus}
    kraken2-build --db !{params.kraken_db} --clean --threads !{task.cpus}
    fi
    cd "!{params.kraken_db}"
    if [[ -d opts.k2d && -d hash.k2d && -d taxo.k2d ]]; then
        kraken_installed=true
    else
        echo "Something went wrong with kraken install. Missing at least one db file"
        echo "ls !{params.kraken_db}"
        ls !{params.kraken_db}
        exit 1
    fi
    '''
}

process writeSampleDetails {
    // gather metadata from file and return
    output:
        path 'sample_*'

    conda "$projectDir/conda/env-artic"

    shell:
    '''
    #!/usr/bin/env python
    import pandas as pd
    # read in metadata
    df = pd.read_csv("!{params.meta}")
    df = df.dropna(subset=["Seq ID"])
    ps_col = [col for col in df.columns if "primer" in col.lower()][0]
    bc_col = [col for col in df.columns if "barcode" in col.lower()][0]

    # merge in scheme info
    df = df.merge(
        pd.read_csv("!{params.scheme_details}"),
        left_on=ps_col,
        right_on="primer_scheme",
        how="left"
    )

    # write out so it can be gathered sample-by-sample
    for i,r in df.iterrows():
        sample_name = r["Seq ID"]

        with open(f"sample_{sample_name}.data",'w') as out:
            out.write(f"sample_name={sample_name}\\n")
            out.write(f"primer_scheme={r[ps_col]}\\n")
            out.write(f"barcode={r[bc_col]}\\n")
            out.write(f"scheme={r['scheme']}\\n")
            out.write(f"scheme_version={r['scheme_version']}\\n")
            out.write(f"min={r['min']}\\n")
            out.write(f"max={r['max']}\\n")
            out.write(f"scheme_dir={r['scheme_dir']}\\n")

    from pathlib import Path
    print("output:",Path(f"sample_{sample_name}.data",'w').resolve())
    '''
}

process sourceSampleData {
    input:
        path sample_data
    output:
        tuple env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(min), env(max), env(scheme_dir)
    // name=`echo "!{sample_data}" | sed 's/sample_//g' | cut -d. -f1`
    shell:
    '''
    # source useful details about each sample
    # echo Sourcing !{sample_data}...
    # if [[ -f !{sample_data} ]]; then echo !{sample_data} exists; 
    # else echo !{sample_data} was not found; fi
    . !{sample_data}
    # echo sourced vars...
    # echo sample_name: $sample_name
    # echo primer_scheme: $primer_scheme
    # echo barcode: $barcode
    # echo scheme: $scheme
    # echo scheme_version: $scheme_version
    # echo min: $min
    # echo max: $max
    # echo scheme_dir: $scheme_dir
    '''
}

// A test to see that all samples and their required metadata have been located
process viewData {
    debug true
    input:
        tuple env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(min), env(max), env(scheme_dir)
    shell:
    '''echo found: $sample_name $primer_scheme $barcode $scheme $scheme_version $min $max $scheme_dir'''
}

process count_raw_reads {
    publishDir "$params.out/raw_read_counts", pattern: '*.counts.csv', mode: 'copy'

    input:
        tuple env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(min), env(max), env(scheme_dir)
    output:
        tuple env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(min), env(max), env(scheme_dir)
        path('*.counts.csv')
    
    shell:
    '''
    let lines=`zcat !{params.fastqs}/${barcode}/* | wc -l`
    reads=$(( lines/4 ))
    echo "${sample_name},$reads" > ${sample_name}.counts.csv
    '''
}

process guppyplex {
    input:
        tuple env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(min), env(max), env(scheme_dir)
    output:
        tuple path('trimmed.fastq.gz'), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)

    conda "$projectDir/conda/env-artic"

    shell:
    '''
    echo "Running guppyplex"
    # run artic guppyplex (trim reads by size)
    artic guppyplex               \
		--skip-quality-check				\
		--min-length ${min}                    \
		--max-length ${max}                    \
		--directory "!{params.fastqs}/${barcode}" \
		--prefix "${sample_name}"             \
		--output "trimmed.fastq" || \


    # zip it up
    echo "Zipping reads"
    gzip 'trimmed.fastq'

    # empty .command.err if no true errors (since any content in .command.err causes nextflow to error out)
    if [[ `cat .command.err | wc -l` -eq 1 && `grep -q Processing .command.err` ]]; then printf '' > .command.err; fi

    echo "Done"
    '''
}


process kraken {
    cpus 16

    publishDir "$params.out/kraken_trimmed", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple path('trimmed.fastq.gz'), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)
        env kraken_installed
    
    output:
        // tuple path('kraken_trimmed.fastq.gz'), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir) //old
        tuple path("*_barcode*.fastq.gz"), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)

    conda "$projectDir/conda/env-kraken2"

    shell:
    '''
    # set thread count
    # if [[ -n $SLURM_CPUS_ON_NODE ]]; then threads=$SLURM_CPUS_ON_NODE; else threads=`nproc`; fi

    echo "running kraken2"
    kraken2 \
        --db "!{params.kraken_db}" \
        --threads !{task.cpus} \
        --unclassified-out "${sample_name}_${barcode}.fastq" \
        'trimmed.fastq.gz'
    
    echo "zipping reads"
    gzip "${sample_name}_${barcode}.fastq"

    echo "checking .command.err"
    [[ $? -eq 0 && `cat .command.err | wc -l` -eq 4 ]] && mv .command.err .kraken.log && touch .command.err

    echo "looks complete"
    '''
}

process porechop {
    cpus 16
    publishDir "$params.out/porechop_kraken_trimmed", mode: 'copy'
    // TODO: set_permissions

    input:
        tuple path('kraken_trimmed.fastq.gz'), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)
    output:
        // tuple path('porechop_kraken_trimmed.fastq.gz'), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir) //old
        tuple path("*_barcode*.fastq.gz"), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)

    conda "$projectDir/conda/env-porechop"

    shell:
    '''
    # set thread count
    # if [[ -n ${SLURM_CPUS_ON_NODE:-} ]]; then threads=${SLURM_CPUS_ON_NODE}; else threads=`nproc`; fi

    porechop \
        -i 'kraken_trimmed.fastq.gz' \
        -o "${sample_name}_${barcode}.fastq.gz" \
        -t !{task.cpus} \
        --check_reads 5000 \
        --min_split_read_size 300 \
            > /dev/null 2>&1
    '''
        // #    | grep -v $'\r' > "${porechop_log_dir}/${run_name}.log"
}

process count_final_reads {
    // publishDir "$params.out/read_counts", pattern: '*.counts.csv', mode: 'copy'

    input:
        tuple path(final_trimmed_fastq), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)
    output:
        tuple path(final_trimmed_fastq), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)
        path('*.counts.csv')

    shell:
    '''
    let lines=`zcat !{final_trimmed_fastq} | wc -l`
    reads=$(( lines/4 ))
    echo "${sample_name},$reads" > ${sample_name}.counts.csv
    '''
}

/*
// // Files produced by minion:
// we keep these
*.consensus.fasta   *.coverage_mask.txt *.primersitereport.txt *.pass.vcf.gz   *.merged.vcf
*.minion.log.txt    *.pass.vcf.gz.tbi   *.primertrimmed.rg.sorted.bam.bai   *.primertrimmed.rg.sorted.bam

// we remove these
*.preconsensus.fasta    *.nCoV-2019_2.vcf   *.nCoV-2019_2.hdf   *.nCoV-2019_1.vcf   *.nCoV-2019_1.hdf   
*.trimmed.rg.sorted.bam.bai *.trimmed.rg.sorted.bam *.sorted.bam.bai    *.sorted.bam    *.primers.vcf   
*.fail.vcf  *.muscle.out.fasta  *.muscle.in.fasta   *.coverage_mask.txt.nCoV-2019_2.depths  
*.coverage_mask.txt.nCoV-2019_1.depths  *.alignreport.txt   *.alignreport.er
*/

process minion {
    // set_permissions
    publishDir "$params.out/samples", pattern: "*.{cons,primersite,pass,merged,primertrimmed,minion}*.*", mode: 'copy'
    publishDir "$params.out/samples", pattern: "*.coverage_mask.txt", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 2
    cpus 16

    input:
        tuple path("final_trimmed.fastq.gz"), env(sample_name), env(primer_scheme), env(barcode), env(scheme), env(scheme_version), env(scheme_dir)
    output:
        // tuple path("*.consensus.fasta"), optional: true, emit: fasta
        tuple path("*.consensus.fasta"), env(sample_name), env(barcode), env(minion_complete)
        path('*')

    conda "$projectDir/conda/env-artic"

    shell:
    '''
    # set thread count
    # if [[ -n ${SLURM_CPUS_ON_NODE:-} ]]; then threads=${SLURM_CPUS_ON_NODE}; else threads=`nproc`; fi

    # run up to 2 times
    minion_complete=false; re=""
	for run in 1 2; do
		if [[ $run = 2 ]]; then re="re"; else re=""; fi
		if [[ ${minion_complete} = false ]]; then
        ( artic minion \
            --medaka \
            --medaka-model r941_prom_high_g360 \
            --no-longshot \
            --normalise 400 \
            --threads !{task.cpus} \
            --scheme-directory ${scheme_dir} \
            --scheme-version ${scheme_version} \
            --read-file "final_trimmed.fastq.gz" \
            "${scheme}" \
            "${sample_name}"
        ) && minion_complete=true || minion_complete=false
        fi
    done

    [[ ${minion_complete} = false ]] && printf '' > ${sample_name}.consensus.fasta
    '''
}

process homopolish {
    cpus 12
    // set_permissions
    publishDir "$params.out/samples", pattern: '*homopolished.fasta', mode: 'copy'
    publishDir "$params.out/samples", pattern: '*final.consensus.fasta', mode: 'copy'

    input:
        tuple path(consensus_fasta), env(sample_name), env(barcode), env(minion_complete)
    output:
        path("*.final.consensus.fasta")
        tuple env(sample_name), env(barcode), env(minion_complete)
        path '*homopolished.fasta'
    // output: homopolish.fasta if created, else consensus.fasta

    conda "$projectDir/conda/env-homopolish"

    shell:
    '''
    # set thread count
    # if [[ -n ${SLURM_CPUS_ON_NODE:-} ]]; then threads=${SLURM_CPUS_ON_NODE}; else threads=`nproc`; fi

    if [[ ${minion_complete} = true ]]; then
        ( homopolish polish \
            -a "!{consensus_fasta}" \
            -t !{task.cpus} \
            -m R9.4.pkl \
            -o . \
            --local_DB_path "!{projectDir}/homopolish_db.fasta"
        ) && homopolish_successful=true || homopolish_successful=false
    else
        touch "${sample_name}_homopolished.fasta"
    fi

    [[ $? -eq 0 ]] && mv .command.err .homopolish.log && touch .command.err

    # choose best consensus to final.consensus.fasta (prioritize homopolish fasta if exists)
    if [[ -f ${sample_name}_homopolished.fasta ]]; then
        echo using homopolished fasta
        best_fasta="${sample_name}_homopolished.fasta"
    elif [[ -f "!{consensus_fasta}" ]]; then
        echo "using artic's fasta"
        best_fasta="!{consensus_fasta}"
    else
        # write empty fasta if none exist
        printf '' > "${sample_name}.final.consensus.fasta"
        return 0
    fi

    # simplify header while copying best fasta to final fasta
    echo ">${sample_name}" > "${sample_name}.final.consensus.fasta"
    sed '1d' "${best_fasta}" >> "${sample_name}.final.consensus.fasta"
    '''
}

process combine_filtered_read_counts {
    // publishDir "$params.out/overall", pattern: 'read_counts*.csv', mode: 'copy', overwrite: true

    input:
        path('*.counts.csv')
    output:
        path('read_counts*.csv')
    
    shell:
    '''
    cat *.counts.csv > read_counts.csv
    '''
}

process combine_raw_read_counts {
    publishDir "$params.out/overall", pattern: 'read_counts*.csv', mode: 'copy'

    input:
        path('*.counts.csv')
    output:
        path('read_counts*.csv')
        path('read_counts*.csv')
    
    shell:
    '''
    cat !{params.out}/raw_read_counts/* > read_counts_raw.csv
    '''
}

process combine_seqs {
    publishDir "$params.out/overall", pattern: 'consensus_*.fasta', mode: 'copy'

    input:
        path '*.fasta'
    output:
        path 'consensus_*.fasta'
        path 'consensus_*.fasta'

    shell:
    '''
    cat *.fasta > "consensus_!{params.plate}.fasta"
    '''
}

process pangolin {
    cpus 16
    publishDir "$params.out/overall", pattern: 'lineage_report-*.csv', mode: 'copy'

    input:
        path fasta
    output:
        path "lineage_report-${params.plate}.csv"

    conda "$projectDir/conda/env-pangolin"

    shell:
    '''
    # update pangolin if needed
    echo updating pangolin if needed    
    !{projectDir}/check4updates.sh pangolin
    [[ $? -eq 0 ]] && mv .command.err pangolin.update.log && touch .command.err

    # set thread count
    # if [[ -n ${SLURM_CPUS_ON_NODE:-} ]]; then threads=${SLURM_CPUS_ON_NODE}; else threads=`nproc`; fi

    echo running pangolin
    pangolin \
        !{fasta} \
        --outfile="lineage_report-!{params.plate}.csv" \
        -t !{task.cpus}

    [[ $? -eq 0 ]] && mv .command.err .pangolin.log && touch .command.err
    '''
}

process nextclade {
    cpus 16
    publishDir "$params.out/overall", pattern: 'nextclade*.csv', mode: 'copy'
    publishDir "$params.out", pattern: 'nextclade', mode: 'copy'

    input:
        path fasta
    output:
        path "nextclade_results-${params.plate}.csv", emit: 'csv'
        path 'nextclade', emit: 'directory'
        env nextclade_version, emit: 'version'

    conda "$projectDir/conda/env-nextclade"

    shell:
    '''
    # update nextclade if needed
    echo updating nextclade if needed
    !{projectDir}/check4updates.sh nextclade
    [[ $? -eq 0 ]] && mv .command.err nextclade.update.log && touch .command.err

    # set thread count
    # if [[ -n ${SLURM_CPUS_ON_NODE:-} ]]; then threads=${SLURM_CPUS_ON_NODE}; else threads=`nproc`; fi

    echo running nextclade
	nextclade run \
		!{fasta} \
		--jobs !{task.cpus} \
		--input-dataset="!{projectDir}/nextclade_data" \
		--output-csv="nextclade_results-!{params.plate}.csv" \
		--output-all "nextclade" \
        --output-selection=fasta,json,ndjson,csv,tree,translations,insertions,errors

    nextclade_version="`nextclade -V | sed 's/nextclade //g'`"
    echo version: $nextclade_version
    '''
}

process create_spreadsheet {
    publishDir "$params.out/overall", pattern: 'Sequencing-report*.csv', mode: 'copy'

    input:
        path filtered_read_counts_csv
        path raw_read_counts_csv
        path pangolin_csv
        path nextclade_csv
        val nextclade_version
    output:
        path 'Sequencing-report*.csv'

    shell:
    '''
    echo running createSpreadsheet.py
    !{projectDir}/createSpreadsheet.py \
        "!{params.meta}" \
        "!{pangolin_csv}" \
        "!{nextclade_csv}" \
        !{nextclade_version} \
        "!{filtered_read_counts_csv}" \
        "!{raw_read_counts_csv}" \
        !{params.plate} \
        . \
        "!{params.reportMap}"
    echo spreadsheet complete
   '''
}

// process checkTest {
//     cpus 1
//     debug true

//     conda "$projectDir/conda/env-pangolin"

//     shell:
//     '''
//     # update pangolin if needed
//     !{projectDir}/check4updates.sh pangolin

//     echo pangolin help
//     pangolin -h

//     '''
// }

workflow {
    // checkTest()
    // Part 1: Artic pipeline
    // ensure kraken db exists, else download
    kraken_complete = check4kraken()
    // Step 1
    // Prep: emit sample metadata to be used individually for each sample
    ch_sample = writeSampleDetails() | flatten
    ch_sample_data = sourceSampleData(ch_sample)

    // Main Pipeline: run individual tasks of artic pipeline
    // viewData(ch_sample_data) // not needed - just to check on things

    // Step 2 count raw reads
    count_raw_reads(ch_sample_data)

    // Step 3
    // artic guppyplex
    guppyplex(count_raw_reads.out[0])

    // Step 4
    kraken(guppyplex.out,kraken_complete)

    // Step 5
    porechop(kraken.out)
    final_trimmed_reads = porechop.out

    // Step 6
    count_final_reads(final_trimmed_reads)

    // Step 7
    minion(count_final_reads.out[0])

    // Step 8
    homopolish(minion.out[0])

    // WAIT FOR ALL FASTAS (for everything above to complete)
    fastas = homopolish.out[0].collect()
    // fastas.view()

    // Step 9
    combine_raw_read_counts(count_raw_reads.out[1].collect())
    combine_filtered_read_counts(count_final_reads.out[1].collect())
    combine_seqs(fastas)

    // Step 10
    pangolin(combine_seqs.out[0])
    nextclade(combine_seqs.out[1])

    // // Step 8a
    // // trigger C-WAP (only ww?) - trigger here and forget?
    // // Step 8b
    // // trigger freyja (only ww?)

    // Step 11
    // Wait on other processes to complete, then create spreadsheet
    create_spreadsheet(combine_filtered_read_counts.out,combine_raw_read_counts.out[1],pangolin.out,nextclade.out.csv,nextclade.out.version)

}


// Lines to remove:
// echo ls:
// ls ...
// debug true
// echo ...
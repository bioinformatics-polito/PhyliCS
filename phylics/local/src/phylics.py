#!/usr/bin/env python

from check_funcs import *
from constants import *
import os
import sys
import glob
import argparse
import multiprocessing

def run_cmd(cmd):
    #print(cmd)
    os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PhyliCS: multi-sample Phylogenetic analysis of Single Cell CNV profiles.")
    run_grp = parser.add_argument_group("execution modes (REQUIRED)", "One of the following execution modes MUST be specified")
    run_cmds = run_grp.add_mutually_exclusive_group(required=True)

    run_cmds.add_argument("--run", action='store_true', help='Run the pipeline. USAGE: phylics --run [--run_cnvs | --run_single | --run_multiple] --input_dirs sample:beds_dir [sample:beds_dir ...] --genome hg19 --binning variable_500000_101_bwa  --meta_format {json,xml} --output_path out_path [--output_prefix out_prefix] [--verbose]') 
    #run_cmds.add_argument("--run_cnv", action='store_true', help='Run only the CNV calling stage. USAGE: phylics.py --input_dirs beds_dir [beds_dir ...] --genome hg19 --binning variable_500000_101_bwa [--init_ginkgo]')
    #run_cmds.add_argument("--run_single", action='store_true', help='Run only the single-sample analysis stage. USAGE: phylics.py --run_single --input_dirs sample:input_dir [sample:input_dir ...]  --meta_format {json,xml} --outdir outdir')
    #run_cmds.add_argument("--run_multiple", action='store_true', help='Run only the  multiple-sample analysis stage. USAGE: phylics.py --run_multiple --input_dirs sample:input_dir [sample:input_dir ...] --meta_format {json,xml} --outdir outdir')
    run_cmds.add_argument("--run_10x_preproc", action='store_true', help='Run 10x data pre-processing module. USAGE: phylics --run_10x_preproc --input_dirs sample_name:10x_out_path --output_path out_oath [--output_prefix out_prefix] [--verbose]. Only single-sample execution is available for this option: only the first input directory is considered, even if more than one has been declared.')
    run_cmds.add_argument("--run_cell_filtering", action='store_true', help='Run the cell filtering module. USAGE: phylics --run_cell_filtering --input_dirs sample_name:input_path --intervals [v1-v2 [v1-v2 ...] --values [v [v ...]] --output_path out_path [--output_prefix out_prefix] [--verbose]. Only single-sample execution is available for this option: only the first input directory is considered, even if more than one has been declared. NOTE that at least one  of the two options, "--intervals" and "--values", must contain values to make this command effective.')

    run_grp1 = parser.add_argument_group("Single-stage execution options", "The following options may be specified to run only one of the pipeline stages when executing the --run mode.")
    run_single_modules_opts = run_grp1.add_mutually_exclusive_group()
    run_single_modules_opts.add_argument("--run_cnvs", action='store_true', help='Run only the CNV calling stage. USAGE: phylics --run --run_cnvs --input_dirs beds_dir [beds_dir ...] --genome hg19 --binning variable_500000_101_bwa [--init_ginkgo] [--verbose]')
    run_single_modules_opts.add_argument("--run_single", action='store_true', help='Run only the single-sample analysis stage. USAGE: phylics --run --run_single --input_dirs sample:input_dir [sample:input_dir ...]  --meta_format {json,xml} --output_path out_path [--output_prefix out_prefix] [--verbose]')
    run_single_modules_opts.add_argument("--run_multiple", action='store_true', help='Run only the  multiple-sample analysis stage. USAGE: phylics --run --run_multiple --input_dirs sample:input_dir [sample:input_dir ...] --meta_format {json,xml} --output_path out_path [--output_prefix out_prefix] [--verbose]')

    run_opts = parser.add_argument_group(title = "execution options")

    #run_opts.add_argument("--sample", required='--run_single' in sys.argv or '--run_cell_filtering' in sys.argv, metavar='sample_name:SegCopy', action='store', help='Sample name and the path to the corresponding cnvs file, separated by ":". Sample name and cnvs file path cannot contain ":"', type=str)
    #run_opts.add_argument("--samples", required=True, metavar='sample_name', help='Sample name list', action='store', nargs='+', type=str)

    run_opts.add_argument("--input_dirs", required=True, metavar='sample_name:input_dir', action=check_valid_samples(), help='Pairs made of sample name and directory path, separated by ":", for each input sample. Sample name and input directory path cannot contain ":"', nargs='+', type=str)

    #run_opts.add_argument("--cnv_files", required='--run_single' in sys.argv or '--run_multiple' in sys.argv or '--run_cell_filtering' in sys.argv, metavar='sample_name:SegCopy', action=check_valid_samples(), help='Sample names and the path to the corresponding cnvs files, separated by ":". Sample names and cnvs file paths cannot contain ":"', nargs='+', type=str)

    #run_opts.add_argument("--other_cnvs", required='--run'  in sys.argv or '--run_multiple' in sys.argv,   metavar='SegCopy', action='store', help='Other CNV filepaths listed in the same order the corresponding sample names have been specified. e.g. --other_samples sample2 sample3 --other_cnvs SegCopy2 SegCopy3', nargs='+', type=str)

    #run_opts.add_argument("--results", required='--run_single' in sys.argv in sys.argv or '--run_cell_filtering' in sys.argv,  metavar='results.txt', action='store', help='Results.txt file path', type=str)

    run_opts.add_argument("--method", choices=['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'], action='store', default='complete', help='Clustering method (default = complete)', type=str)
    run_opts.add_argument("--metric", choices=['euclidean', 'cityblock', 'sqeuclidean', 'cosine', 'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra'], default='euclidean', action='store', help='Distance metric', type=str)
    run_opts.add_argument("--meta_format", required='--run'  in sys.argv, choices=['json', 'xml'], action='store', help='Metadata file format.', type=str)

    #run_opts.add_argument("--outdirs", required='--run_single' in sys.argv or '--run_cell_filtering' in sys.argv, metavar=':sample:outdir', action=check_valid_samples(), help='Sample names and the path to the corresponding output directory, separated by ":". Sample names and output directory paths cannot contain ":".', nargs='+', type=str)

    run_opts.add_argument("--output_path", required='--run_cnvs' not in sys.argv, metavar='out_path', action='store', help='Path to the location where the output directories for the different analysis stages must be created.',  type=str)
    run_opts.add_argument("--output_prefix", metavar='out_prefix', action='store', help='A string to be pre-pended to the names automatically generated for the output directories',  type=str)
    run_opts.add_argument("--tasks", metavar='N', action='store', help='Maximum number of tasks to be run in parallel. It allows to execute single-sample analyses in parallel and to parallelize the permutation test execution for the heterogeneity score computation.',
                                               type=int)
    run_opts.add_argument("--seed", metavar='N', 
                help='Seed to initialize the pseudo-random generator used to perform the permutation test.', type=int)
    
    run_opts.add_argument("--n_permutations", metavar='N', action='store', 
                        help='Number of permutations to execute the permutation test for the heterogeneity score computation.',
                        type=int)
    '''
    run_opts.add_argument("--tsne_iterations", metavar='N', action='store', 
        help='Number of iterations for tSNE computation.',
                        type=int)
    run_opts.add_argument("--tsne_perplexity", metavar='N', action='store', 
                help='Perplexity value for tSNE computation.',
                       type=int)
    '''
    run_opts.add_argument("--reclust",metavar='N', action=check_valid_N_clust(), 
                        help='If this option is specified, only the clustering part is executed with the specified number of clusters (see --n_clusters option), unless --reinit option is specified (see below).',
                                                type=int)
    run_opts.add_argument("--reinit", action='store_true',
                help='This option has effect only if combined with the --reclust option. It allows to recompute the entire analysis and then recluster with the specified number of clusters.')
    run_opts.add_argument('--intervals', action=valid_interval(), metavar='v1-v2', help='List of of mean ploidy intervals: cells which mean ploidies are in the specified ranges are filtered out', nargs='+')
    run_opts.add_argument('--values',  action='store', metavar='v', help='List of of mean ploidy values: cells which mean ploidies are equal to the specified ones are filtered out', nargs='+', type=int)

    #run_opts.add_argument('--beds_dir', required='--run_cnv' in sys.argv, action='store', metavar='path/to/beds/dir', help='Path to the directory where compressed .bed files are located', type=str)

    run_opts.add_argument('--genome', required='--run_cnvs' in sys.argv and '--run' in sys.argv, action='store', metavar='genome', help='Directory name for ROOT_DIR/genomes/${chosen_genome}', type=str)
    run_opts.add_argument('--binning',  required='--run_cnvs' in sys.argv and '--run' in sys.argv, action='store', metavar='variable_500000_101_bwa', help='A complex value made of the concatenation of - type: variable or fixed (bins. Variable refers to amount of mappable genome, recommended); - size: bin size; - read-length: mapped reads length; - aligner: bowtie or bwa. The read-length and aligner refer to the simulations of re-mapping reads of that length with that aligner on the whole genome. This is used to calculate bins of "mappable" (i.e. variable) genome.The resulting value is the name of a file under ginkgo/genomes/$choosen_genome/original/ with the bin coordinates.', type=str)
    #run_opts.add_argument('--init_ginkgo', action='store_true', help='Clean the directory and start from scratch the whole analysis')
    
    #run_opts.add_argument("--barcodes_csv", required = '--run_10x_preproc' in sys.argv, action='store', metavar='per_cell_summary_metrics.csv', help='Cell Ranger DNA output file containing cell barcodes.', type=str)

    #run_opts.add_argument("--input_bam", required = '--run_10x_preproc' in sys.argv, action='store', metavar='possorted_bam.bam', help='Cell Ranger DNA output file containing reads alignments.' type=str)
    
    run_opts.add_argument("--verbose", action='store_true',  help='Verbose execution.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    '''
    elif '--run' not in sys.argv and '--run_10x_preproc' not in sys.argv and '--run_cell_filtering' not in sys.argv and '--help' not in sys.argv and '-h' not in sys.argv:
        print("ERROR: one of the following execution modes must be specified: --run, --run_10x_preproc, --run_cell_filtering")
        parser.print_help()
        sys.exit(1)
    '''
    args = parser.parse_args()
 

    if args.run:
        #print("Complete execution")
        #stages: cnv calling, single-sample analysis
        #print('ArgumentError: --input_dirs requires a list of items of the same size of --samples')
        #sys.exit(1)
        run_cnvs = True
        run_single = True
        run_multiple = True

        if args.run_cnvs:
            run_single = False
            run_multiple = False
        elif args.run_single:
            run_cnvs = False
            run_multiple = False
        elif args.run_multiple:
            run_cnvs = False
            run_single = False

        input_dirs = {}
        for sample_input in args.input_dirs:
            sample = sample_input.split(':')[0]
            inputdir = sample_input.split(':')[1]
            if os.path.isdir(inputdir):
                input_dirs[sample] = inputdir
            else:
                print("FilePathError: {} is not a valid path".format(inputdir))
                sys.exit(1)
        
        out_path = ""
        out_prefix =  ""
        
    
        samples = (input_dirs.keys())
        #if sorted(args.samples) != sorted(sample_inputs.keys()):
            #print('ArgumentError: --input_dirs and --samples must specify the same set of samples')
        if run_cnvs:
            #-------------------------------------#
            #           cnvs calling              #
            #-------------------------------------#

            cnvs_cmds = []
            for sample in samples:
                input_dir = os.path.join(os.getcwd(), input_dirs[sample])
        
                ginkgo_path = os.path.join(BIN_DIR,"ginkgo", "cli", "ginkgo.sh")
                cmd = "{} --input {} --genome {} --binning {} --init | grep -v \"pipe\"" .format(ginkgo_path, input_dir, args.genome, args.binning)
                #if args.init_ginkgo:
                #    cmd = cmd + " --init"
                cnvs_cmds.append(cmd)

            if args.tasks:
                njobs = args.tasks
                nsamples = len(samples)

                for i in range(0, nsamples, njobs):
                    jobs = []
                    for j in range(0, njobs):
                        if i < nsamples:
                            cmd = cnvs_cmds[i]
                            #print("Sample {} processing starting".format(samples[i])
                            process = multiprocessing.Process(target=run_cmd,args=(cmd,))
                            jobs.append(process)
                            i = i + 1
                    # Start the processes 
                    for job in jobs:
                        job.start()
                    # Ensure all of the processes have finished
                    for job in jobs:
                        job.join()
                    #print("{} jobs done.".format(njobs))
            else:
                for cmd in cnvs_cmds:
                    if args.verbose:
                        os.system(cmd)
                        #print(cmd)
                    else:
                        cmd = cmd + "  &> ginkgo.log"
                        #print(cmd)
                        os.system(cmd)

        if run_single:
            #-------------------------------------#
            #           single-sample             #
            #              analysis               #
            #-------------------------------------#
            out_path = args.output_path
            if args.output_prefix:
                out_prefix = args.output_prefix


            cmds = []
            for sample in samples:
                #sample = sample_cnvs.split(":")[0]
                cnvs_file = os.path.join(input_dirs[sample], "SegCopy")
                if os.path.exists(cnvs_file) == False:
                    print("FilePathError: {} does not  exist".format(cnvs_file))
                    sys.exit(1)
                results_file = os.path.join(input_dirs[sample], "results.txt")
                if os.path.exists(results_file) == False:
                    print("FilePathError: {} does not  exist".format(results_file))
                    sys.exit(1)
                
                out_name = sample + "_post_CNV"
                if out_prefix != "":
                    out_name = out_prefix + "_" + out_name
                out_dir= os.path.join(out_path, out_name)
                if os.path.exists(out_dir) == False:
                    try:
                        os.mkdir(out_dir)
                    except OSError:
                        print ("Creation of the directory {} failed".format(out_dir))
                        sys.exit(1)
                    else:
                        print ("Successfully created the directory {} ".format(out_dir))

                tool = "single_sample_post_analysis"
                cmd = "{} {} {} {} {} {} {} {}".format(tool, sample, cnvs_file, results_file, args.method, args.metric, args.meta_format, out_dir)
                '''
                if args.tsne_iterations:
                    cmd = "{} --tsne_iterations {}".format(cmd, args.tsne_iterations)
                if args.tsne_perplexity:
                    cmd = "{} --tsne_perplexity {}".format(cmd, args.tsne_perplexity)
                '''
                if args.reclust:
                    cmd = "{} --reclust {}".format(cmd, args.reclust)
                if args.reinit:
                    cmd = "{} --reinit".format(cmd)
                if args.verbose:
                    cmd = "{} --verbose".format(cmd)
                cmds.append(cmd)
        
            if args.tasks:
                njobs = args.tasks
                nsamples = len(samples)

                for i in range(0, nsamples, njobs):
                    jobs = []
                    for j in range(0, njobs):
                        if i < nsamples:
                            cmd = cmds[i]
                            #print("Sample {} processing starting".format(samples[i])
                            process = multiprocessing.Process(target=run_cmd,args=(cmd,))
                            jobs.append(process)
                            i = i + 1
                    # Start the processes 
                    for job in jobs:
                        job.start()
                    # Ensure all of the processes have finished
                    for job in jobs:
                        job.join()
                    #print("{} jobs done.".format(njobs))
            else:
                for cmd in cmds:
                    #print(cmd)
                    os.system(cmd)
        if run_multiple:
            #-------------------------------------#
            #           multiple-sample           #
            #              analysis               #
            #-------------------------------------#
            if out_path == "":
                out_path = args.output_path

            if args.output_prefix and out_prefix == "":
                out_prefix = args.output_prefix
            
            out_name = ""
            for sample in samples:
                if out_name == "":
                    out_name = sample
                else:
                    out_name = out_name + "_" + sample
            
            if out_prefix != "":
                out_name = out_prefix + "_" + out_name
            out_name = out_name + "_postCNV"
            out_dir = os.path.join(out_path, out_name)
            if os.path.exists(out_dir) == False:
                try:
                    os.mkdir(out_dir)
                except OSError:
                    print ("Creation of the directory {} failed".format(out_dir))
                    sys.exit(1)
                else:
                    print ("Successfully created the directory {} ".format(out_dir))

            cnvs_files = []
            for sample in samples:
                cnvs_file = os.path.join(input_dirs[sample], "SegCopy")
                if os.path.exists(cnvs_file) == False:
                    print("FilePathError: {} does not  exist".format(cnvs_file))
                    sys.exit(1)
                sample_cnvs = sample + ":" + cnvs_file
                cnvs_files.append(sample_cnvs)  

            sample_cnvs_files  = ' '.join(cnvs_files)
                        
            tool =  "multi_sample_post_analysis"
            cmd = "{} {} {} {} {} {}".format(tool, sample_cnvs_files, args.method, args.metric, args.meta_format, out_dir)
            if args.tasks:
                cmd = "{} --n_jobs {}".format(cmd, args.tasks)
            if args.n_permutations:
                cmd = "{} --n_permutations {}".format(cmd, args.n_permutations)
            '''
            if args.tsne_iterations:
                cmd = "{} --tsne_iterations {}".format(cmd, args.tsne_iterations)
            if args.tsne_perplexity:
                cmd = "{} --tsne_perplexity {}".format(cmd, args.tsne_perplexity)
            '''
            if args.reclust:
                cmd = "{} --reclust {}".format(cmd, args.reclust)
            if args.reinit:
                cmd = "{} --reinit".format(cmd)
            if args.verbose:
                cmd = "{} --verbose".format(cmd)
            #print(cmd)
            os.system(cmd)

    elif args.run_10x_preproc:
        print("10x data preprocessing execution")
        """
           {params.tool} --barcodes-csv {input.csv} --forbidden-tags XA,SA --min-mapq 30 -o {output} {input.bam}  &> {log}

            mv {input}/noise.bam {params.noisedir}/{params.newnoisebam}
                
            samtools view  -u {input} | bamToBed -i - > {output}
          
            gzip {input}
                
        """
        if os.path.exists(args.output_path) == False:
            print("FilePathError: {} does not  exist".format(args.output_path))
            sys.exit(1)


        sample_input = args.input_dirs[0]
        sample = sample_input.split(":")[0]
        input_dir = sample_input.split(":")[1]

        input_bam = os.path.join(input_dir, "possorted_bam.bam")
        if os.path.exists(input_bam) == False:
            print("FilePathError: {} does not  exist".format(input_bam))
            sys.exit(1)
        barcodes_csv = os.path.join(input_dir, "per_cell_summary_metrics.csv")
        if os.path.exists(barcodes_csv) == False:
            print("FilePathError: {} does not  exist".format(barcodes_csv))
            sys.exit(1)


        out_name = sample + "_demux"
        out_path = args.output_path

        if args.output_prefix :
            out_name = args.output_prefix + "_" + out_name

        out_dir = os.path.join(out_path, out_name)

        if os.path.exists(out_dir) == False:
                try:
                    os.mkdir(out_dir)
                except OSError:
                    print ("Creation of the directory {} failed".format(out_dir))
                    sys.exit(1)
                else:
                    print ("Successfully created the directory {} ".format(out_dir))

        demux = "sctools_demultiplex"
        log_file = os.path.join(out_dir, "demux.log")
    
        
        if args.verbose:
            os.system("{tool} --barcodes-csv {csv} --forbidden-tags XA,SA --min-mapq 30 -o {output} {bam}".format(tool=demux, csv=barcodes_csv, output=out_dir, bam=input_bam))
        else:
            os.system("{tool} --barcodes-csv {csv} --forbidden-tags XA,SA --min-mapq 30 -o {output} {bam} &> {log}".format(tool=demux, csv=barcodes_csv, output=out_dir, bam=input_bam, log=log_file))
        
        # mv {input}/noise.bam {params.noisedir}/{params.newnoisebam}
        noise_bam = os.path.join(out_dir, "noise.bam")
        noise_rename = os.path.join(out_dir, "noise.bam.noise")
        os.rename(noise_bam, noise_rename)
        
        # samtools view  -u {input} | bamToBed -i - > {output}
        bams = []
        for file in glob.glob(os.path.join(out_dir, "*.bam")):
            bams.append(file)

        
        cmds = []
        for bam in bams:
            bed = os.path.splitext(bam)[0] + ".bed"
            gzip = bed + ".gz"
            cmds.append("samtools view  -u {bam_file} | bamToBed -i - > {bed_file} | gzip > {gzipped}".format(bam_file=bam, bed_file=bed, gzipped=gzip))
        
        if args.tasks:
            njobs = args.tasks
            nbams = len(bams)
                
            for i in range(0, nbams, njobs):
                jobs = []
                for j in range(0, njobs):
                    if i < nbams:
                        cmd = cmds[i]
                        process = multiprocessing.Process(target=run_cmd,args=(cmd,))
                        jobs.append(process)
                        i = i + 1
                # Start the processes 
                for job in jobs:
                    job.start()
                # Ensure all of the processes have finished
                for job in jobs:
                    job.join()
                #print("{} jobs done.".format(njobs))
        else:
            for cmd in cmds:
                if args.verbose:
                    os.system(cmd)
                    #print(cmd)
                else:
                    cmd = cmd + "  &> samtools.log"
                    #print(cmd)
                    os.system(cmd)

        
    elif args.run_cell_filtering:
        print("Cell filtering execution")
        """
        usage: valid_cells.py [-h]
                          sample_name results.txt SegCopy interval/value
                          [interval/value ...] out_dir
        """
        if os.path.exists(args.output_path) == False:
            print("FilePathError: {} does not  exist".format(args.output_path))
            sys.exit(1)

        
        sample_input = args.input_dirs[0]
        sample = sample_input.split(":")[0]
        input_dir = sample_input.split(":")[1]

        cnvs = os.path.join(input_dir, "SegCopy")
        if os.path.exists(cnvs) == False:
            print("FilePathError: {} does not  exist".format(cnvs))
            sys.exit(1)
        results = os.path.join(input_dir, "results.txt")
        if os.path.exists(results) == False:
            print("FilePathError: {} does not  exist".format(results))
            sys.exit(1)     

        if args.values == None and args.intervals == None:
            print("phylics: error: expected at least one of the following arguments: --values, --intervals")
            sys.exit(3)

        out_name = sample + "_filtered"
        out_path = args.output_path
        
        if args.output_prefix :
            out_name = args.output_prefix + "_" + out_name

        out_dir = os.path.join(out_path, out_name)

        if os.path.exists(out_dir) == False:
                try:
                    os.mkdir(out_dir)
                except OSError:
                    print ("Creation of the directory {} failed".format(out_dir))
                    sys.exit(1)
                else:
                    print ("Successfully created the directory {} ".format(out_dir))


        intervals = ""
        values = ""

        if args.intervals:
            if (len(args.intervals) != 0):
                intervals = " ".join(args.intervals)
        if args.values:
            if (len(args.values) !=0):
                values = " ".join(args.values)
    
        outliers = intervals + " " + values

        tool = os.path.join("valid_cells")
        cmd = "{} {} {} {} {} {}".format(tool, sample, results, cnvs, outliers, out_dir)
        if args.verbose:
            cmd = "{} --verbose".format(cmd)
    
        #print(cmd)
        os.system(cmd)

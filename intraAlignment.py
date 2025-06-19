# -*- coding: utf-8 -*-
"""
Project         : singleMicrobiome
File            : intraAlignment.py
Creator         : Ping
Create time     : 2024/08/02
Software        : PyCharm
Introduction    :
"""
import operator, os, sys
import gzip, json, time
import subprocess
import numpy as np
import pandas as pd
from colorama import Fore, Style, init
from tqdm import tqdm
from multiprocessing import Pool
import argparse as ap
from functools import wraps

# set working directory to script path
os.chdir(sys.path[0])
# set colorama auto reset
init(autoreset=True)


########################################################
# time_it decorator
########################################################
def time_it(func):
    """
    Decorator to calculate the running time of a function
    :param func:
    :return:
    """

    @wraps(func)
    def decorate(*args, **kwargs):
        start = time.time()
        if func.__name__ != "run":
            print(Fore.LIGHTGREEN_EX + 'Running Function: ' +
                  Fore.LIGHTRED_EX + '%s' % func.__name__)  # .split("_", 1)[-1])
        func(*args, **kwargs)
        end = time.time()
        if func.__name__ == "run":
            print(Fore.LIGHTGREEN_EX + 'Total cost time : %.3f s' % (end - start))
        else:
            print(Fore.LIGHTGREEN_EX + 'Cost time : %.3f s' % (end - start))

    return decorate


########################################################
# SoftwareNotFoundError
########################################################

class intraAlignmentError(Exception):
    """
    Process Initialization Error
    """

    def __init__(self, error_info):
        super().__init__(self)
        self.error_info = error_info
        self.intraAlignment_error = []

    def __str__(self):
        if self.error_info["software"]:
            self.intraAlignment_error.append("%s%s  SoftwareNotFoundError:%s\n%s" % (Fore.LIGHTRED_EX,
                                                                                     Style.BRIGHT,
                                                                                     Style.RESET_ALL,
                                                                                     "\n".join(
                                                                                         self.error_info["software"])))

        if self.error_info["data"]:
            self.intraAlignment_error.append("%s%s  DataNotFoundError:%s\n%s" % (Fore.LIGHTRED_EX,
                                                                                 Style.BRIGHT,
                                                                                 Style.RESET_ALL,
                                                                                 "\n".join(self.error_info["data"])))
        if self.error_info["database"]:
            self.intraAlignment_error.append("%s%s  DatabaseNotFoundError:%s\n%s" % (Fore.LIGHTRED_EX,
                                                                                     Style.BRIGHT,
                                                                                     Style.RESET_ALL,
                                                                                     "\n".join(
                                                                                         self.error_info["database"])))

        return Fore.LIGHTRED_EX + Style.BRIGHT + "Process Initialization Error:\n" + Style.RESET_ALL + "\n".join(
            self.intraAlignment_error)


# check software
def check_software(software_check_command):
    try:
        subprocess.check_output(software_check_command,
                                stderr=subprocess.STDOUT,
                                shell=True)
        return ""
    except subprocess.CalledProcessError:
        return "        %s%-8s%s not installed." % (Fore.LIGHTGREEN_EX,
                                                    software_check_command.split(" ")[0],
                                                    Style.RESET_ALL)


########################################################
# Tree Class, referenced from "https://github.com/jenniferlu717/KrakenTools/blob/master/make_kreport.py"
########################################################

class Tree(object):
    """
    Tree node.
    """

    def __init__(self, tax_id, name, level_rank, level_num, p_tax_id, parent=None, children=None):
        self.tax_id = tax_id
        self.name = name
        self.level_rank = level_rank
        self.level_num = int(level_num)
        self.p_tax_id = p_tax_id
        self.all_reads = 0
        self.lvl_reads = 0
        # Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)


########################################################
# intraAlignment alignment pipeline
########################################################
class intraAlignment:
    def __init__(self,
                 ranger_out_path,
                 align_method,
                 platform,
                 c4_output_path,
                 c4_cdna_r1,
                 reference_database_path,
                 kraken_thread,
                 single_thread,
                 blast_thread,
                 test_mode,
                 clean_mode):

        # set print style
        self.step_style = Fore.LIGHTRED_EX + Style.BRIGHT + "%s"

        # check errors
        self.errors = {"software": [], "database": [], "data": []}
        # set parameters
        self.test_mode = test_mode
        self.clean_mode = clean_mode
        self.ranger_out_path = ranger_out_path
        self.align_method = align_method
        self.reference_database_path = reference_database_path.rstrip("/")

        self.platform = platform
        if platform == "c4":
            self.c4_cdna_r1 = c4_cdna_r1
            self.c4_output_path = c4_output_path
            self.ranger_out_path = c4_output_path.rstrip("/") + "/output"
        # check 10x input data
        elif platform == "10x":
            if not os.path.exists("%s/possorted_genome_bam.bam" % self.ranger_out_path):
                self.errors["data"].append(
                    '        -rp [--rgPath]:       "%s%s%s"\n' % (Fore.LIGHTGREEN_EX,
                                                                  self.ranger_out_path,
                                                                  Style.RESET_ALL) +
                    '        Necessary input file: "%s%s/%spossorted_genome_bam.bam%s" does not exist.' % (
                        Fore.LIGHTGREEN_EX,
                        self.ranger_out_path,
                        Fore.LIGHTRED_EX,
                        Style.RESET_ALL))

        # set output path
        self.output_path = "%s/intraAlignment" % self.ranger_out_path

        # init logger， monitor the pipeline progress
        if os.path.exists("%s/intraAlignment.log" % self.output_path):
            with open("%s/intraAlignment.log" % self.output_path, "r") as f:
                self.logger = json.load(f)
        else:
            self.logger = {}

        if align_method == "kraken2":
            # detect software
            for software in ["samtools view", "kraken2 --version", "bracken"]:
                error = check_software(software)
                if error:
                    self.errors["software"].append(error)
            # detect database
            if not os.path.exists("%s/ktaxonomy.tsv" % self.reference_database_path):
                self.errors["database"].append(
                    '        -db [--dbPath]:              "%s%s%s"\n' % (Fore.LIGHTGREEN_EX,
                                                                         self.reference_database_path,
                                                                         Style.RESET_ALL) +
                    '        kraken2 reference database: ["%s%s/%sktaxonomy.tsv%s"\n'
                    '                                     "%s%s/%s[hash ... taxo].k2d%s"] does not exist.' % (
                        Fore.LIGHTGREEN_EX,
                        self.reference_database_path,
                        Fore.LIGHTRED_EX,
                        Style.RESET_ALL,
                        Fore.LIGHTGREEN_EX,
                        self.reference_database_path,
                        Fore.LIGHTRED_EX,
                        Style.RESET_ALL)
                )

            if self.errors["software"] or self.errors["database"] or self.errors["data"]:
                raise intraAlignmentError(self.errors)
            os.makedirs(self.output_path, exist_ok=True)
            self.output_path = "%s/kraken2" % self.output_path
            os.makedirs(self.output_path, exist_ok=True)
            self.kraken_thread = kraken_thread
            self.single_thread = single_thread

            self.kraken_output_path = "%s/singleCell_kraken_output" % self.output_path
            self.kraken_report_path = "%s/singleCell_kraken_report" % self.output_path
            self.bracken_report_path = "%s/singleCell_bracken_report" % self.output_path
            self.bracken_result_path = "%s/singleCell_bracken_result" % self.output_path
            self.mpa_report_path = "%s/singleCell_mpa_report" % self.output_path
            # create output dirs
            for path in [self.kraken_output_path, self.kraken_report_path,
                         self.bracken_report_path, self.bracken_result_path, self.mpa_report_path]:
                os.makedirs(path, exist_ok=True)
        elif align_method == "blast":
            # detect software
            for software in ["samtools view", "blastn -h", "taxonkit -h"]:
                error = check_software(software)
                if error:
                    self.errors["software"].append(error)
            # detect database
            if sum([os.path.exists("%s/%s" % (self.reference_database_path, x)) for x in
                    ["nt.000.nsq", "names.dmp"]]) != 2:
                self.errors["database"].append(
                    '        -db [--dbPath]:            "%s%s%s"\n' % (Fore.LIGHTGREEN_EX,
                                                                       self.reference_database_path,
                                                                       Style.RESET_ALL) +
                    '        blast reference database: ["%s%s/%snt.###.[nhd ... nsq]%s",\n'
                    '                                   "%s%s/%s[names ... nodes].dmp%s"] does not exists.' % (
                        Fore.LIGHTGREEN_EX,
                        self.reference_database_path,
                        Fore.LIGHTRED_EX,
                        Style.RESET_ALL,
                        Fore.LIGHTGREEN_EX,
                        self.reference_database_path,
                        Fore.LIGHTRED_EX,
                        Style.RESET_ALL)
                )

            if self.errors["software"] or self.errors["database"] or self.errors["data"]:
                raise intraAlignmentError(self.errors)

            os.makedirs(self.output_path, exist_ok=True)
            self.output_path = "%s/blast" % self.output_path
            os.makedirs(self.output_path, exist_ok=True)
            self.blast_thread = blast_thread

    ########################################################
    # write logger
    ########################################################

    def write_logger(func):
        """
        Decorator to write the logger
        :param func:
        :return:
        """

        @wraps(func)
        def decorate(self, *args, **kwargs):
            if func.__name__ not in self.logger:
                func(self, *args, **kwargs)
                self.logger[func.__name__] = "Done"
                with open("%s/intraAlignment/intraAlignment.log" % self.ranger_out_path, "w") as fw:
                    # formatted json file
                    json.dump(self.logger, fw, indent=4)
            else:
                print(Fore.LIGHTRED_EX + "%s has been done, skip this step" % func.__name__)

        return decorate

    ########################################################
    ########################################################
    # c4 platform
    ########################################################
    ########################################################
    @time_it
    def c4_extract_unmapped_reads(self):
        # step1 sorting unmapped fq by sequence labels
        command1 = f"seqkit sort {self.c4_output_path}/01.data/Unmapped.out.mate1 > {self.output_path}/unmapped.sorted.fq"

        # step2. extract seq labels from unmapped fq
        command2 = f"sed -n '1~4s/^@\\([^[:space:]]*\\).*/\\1\\/1/p' {self.output_path}/unmapped.sorted.fq > {self.output_path}/unmapped.labels.txt"

        # step3. extract R1 fastq data by unmapped labels
        if self.c4_cdna_r1.endswith("gz"):
            command3 = (f"zcat {self.c4_cdna_r1} | seqkit grep -f {self.output_path}/unmapped.labels.txt - "
                        f"| seqkit sort - > {self.output_path}/unmapped.r1.fq")
        else:
            command3 = (f"cat {self.c4_cdna_r1} | seqkit grep -f {self.output_path}/unmapped.labels.txt - "
                        f"| seqkit sort - > {self.output_path}/unmapped.r1.fq")

        if not os.path.exists(f"{self.output_path}/unmapped.sorted.fq"):
            print(Fore.LIGHTBLUE_EX + command1 + Style.RESET_ALL)
            subprocess.run(command1, shell=True)
        print(Fore.LIGHTBLUE_EX + command2 + Style.RESET_ALL)
        subprocess.run(command2, shell=True)
        print(Fore.LIGHTBLUE_EX + command3 + Style.RESET_ALL)
        subprocess.run(command3, shell=True)

        # step4. extract unmapped reads with barcode and umi info
        # read output/filtered_feature_bc_matrix/barcodes.tsv.gz
        with gzip.open(f"{self.c4_output_path}/output/filtered_feature_bc_matrix/barcodes.tsv.gz", "rt") as f:
            cell_barcode = set([x.strip() for x in f.readlines()])
            sc_bc = pd.read_csv(f"{self.c4_output_path}/output/singlecell.csv")
            barcode_info = [x.split(";") for x in sc_bc["BARCODE"]]
            barcode_info_filtered = [x for x in barcode_info if x[0] in cell_barcode]
            barcode_dict = {y: x[0] for x in barcode_info_filtered for y in x}
            barcode_len = len(list(barcode_dict.values())[0])

        if self.align_method == "kraken2":
            unmapped_file_name = f"{self.output_path}//possorted_genome_bam_unmapped_noHost.fastq"
        elif self.align_method == "blast":
            unmapped_file_name = f"{self.output_path}//possorted_genome_bam_unmapped_noHost.fasta"
        else:
            raise ValueError(f"Unsupported align method: {self.align_method}")

        with open(f"{self.output_path}/unmapped.sorted.fq") as f_unmapped:
            with open(f"{self.output_path}/unmapped.r1.fq") as f_r1:
                with open(unmapped_file_name, "w") as f_unmapped_out:
                    while True:
                        # read seqs
                        unmapped_lines = [f_unmapped.readline() for _ in range(4)]
                        r1_lines = [f_r1.readline() for _ in range(4)]

                        if not unmapped_lines[0]:
                            break

                        unmapped_label = unmapped_lines[0].split()[0]
                        r1_label = r1_lines[0].split("/")[0]
                        if unmapped_label != r1_label:
                            # delete output/intraAlignment/kraken2/possorted_genome_bam_unmapped_noHost.fastq
                            print("Removing unmapped reconstructed file: possorted_genome_bam_unmapped_noHost.fastq")
                            # os.remove("output/intraAlignment/kraken2/possorted_genome_bam_unmapped_noHost.fastq")
                            raise ValueError(f"Unmapped label {unmapped_label} does not match R1 label {r1_label}")
                        barcode = r1_lines[1].strip()[:barcode_len]
                        umi = r1_lines[1].strip()[barcode_len:]
                        if barcode in barcode_dict:
                            # modify labels with barcode and umi
                            unmapped_lines[0] = "%s|CB:Z:%s|UB:Z:%s" % (unmapped_label, barcode_dict[barcode],
                                                                        umi) + "\n"
                            print(barcode, unmapped_label)
                            if self.align_method == "kraken2":
                                f_unmapped_out.write("".join(unmapped_lines))
                            elif self.align_method == "blast":
                                f_unmapped_out.write(unmapped_lines[0].replace("@", ">") + unmapped_lines[1])

    ########################################################
    ########################################################
    # kraken2 microbial alignment progress
    ########################################################
    ########################################################
    @time_it
    @write_logger
    def kraken2_step1_extract_unmapped_reads(self):
        """
        Extract unmapped reads from the bam file
        :return:
        """

        # check test mode
        if self.test_mode:
            return

        # extract barcode information from bam file
        samtools_command = [
            "samtools view",
            "-@ 24",
            "-f 4",
            "%s/possorted_genome_bam.bam" % self.ranger_out_path
        ]
        awk_command = r'''awk 'match($0,/CB:Z:[ATCG]*-1/,b)&&match($0,/UB:Z:[ATCG]*/,u){print "@"$1"|"b[0]"|"u[0]"\n"$10"\n+"$1"|"b[0]"|"u[0]"\n"$11}' ''' + \
                      '''> %s/possorted_genome_bam_unmapped_noHost.fastq''' % self.output_path

        # Combine the commands into a single string
        extract_command = " ".join(samtools_command) + " | " + awk_command

        subprocess.run(extract_command, shell=True)

    @time_it
    @write_logger
    def kraken2_step2_alignment(self):
        """
        kraken2 microbial alignment progress
        :return:
        """

        # check test mode
        if self.test_mode:
            return

        # check kraken2 alignment result
        if os.path.exists("%s/kraken.report.txt" % self.output_path):
            # print(Fore.LIGHTGREEN_EX + "Kraken2 alignment result already exists, skip this step")
            return
        else:
            subprocess.run(["kraken2",
                            "--db", self.reference_database_path,
                            "--threads", str(self.kraken_thread),
                            "--report", "%s/kraken.report.txt" % self.output_path,
                            "--output", "%s/unmapped_noHost_kraken2.output" % self.output_path,
                            "%s/possorted_genome_bam_unmapped_noHost.fastq" % self.output_path])

    # split kraken2 standard output by barcode
    @time_it
    @write_logger
    def kraken2_step3_split_kraken_output(self):
        """
        Split kraken standard output by barcode
        :return:
        """
        # print(self.step_style % "Step3: Split kraken standard output by barcode")
        if self.test_mode:
            return
        # filter kraken2 output by filtered barcodes list from cellranger output or spaceranger output

        command_1 = "zcat %s/filtered_feature_bc_matrix/barcodes.tsv.gz | grep -Ff - %s/unmapped_noHost_kraken2.output  > %s/unmapped_noHost_kraken2_filtered.output" \
                    % (self.ranger_out_path, self.output_path, self.output_path)
        command_2 = """awk -i inplace '$1 == "C" && $3 != 9606' %s/unmapped_noHost_kraken2_filtered.output""" % self.output_path
        subprocess.run(command_1, shell=True)
        subprocess.run(command_2, shell=True)
        data_dict = {}

        with gzip.open("%s/filtered_feature_bc_matrix/barcodes.tsv.gz" % self.ranger_out_path, "rt") as f:
            barcodes_ranger = [x.strip() for x in f.readlines()]
        total_line_num = \
            subprocess.check_output(
                "wc -l %s/unmapped_noHost_kraken2_filtered.output" % self.output_path,
                shell=True).decode().split(" ")[0]
        line_num = 0

        with open("%s/unmapped_noHost_kraken2_filtered.output" % self.output_path, "r") as f:
            for line in f:
                line_num += 1
                barcode = line.split("CB:Z:")[-1].split("|UB:Z:")[0]
                umi = line.split("UB:Z:")[-1].split("\t")[0]
                if barcode in barcodes_ranger:
                    if barcode not in data_dict.keys():
                        # print("%d / %s %s %s" % (line_num, total_line_num, barcode, umi), end="\r")
                        data_dict[barcode] = {"umi": [umi], "output": line}
                    else:
                        if umi not in data_dict[barcode]["umi"]:
                            # print("%d / %s %s %s" % (line_num, total_line_num, barcode, umi), end="\r")
                            data_dict[barcode]["umi"].append(umi)
                            data_dict[barcode]["output"] += line
        # store intraAlignment_analysis cell kraken output
        for bc in data_dict:
            with open("%s/%s.output" % (self.kraken_output_path, bc), "w") as fa:
                fa.write(data_dict[bc]["output"])

    # make kraken report for intraAlignment_analysis cell
    ########################################################
    # refer to "https://github.com/jenniferlu717/KrakenTools/blob/master/make_kreport.py"
    ########################################################
    def make_kreport_report(self, file):
        """
        Make kraken report for intraAlignment_analysis cell
        :param file:
        :return:
        """
        with open("%s/%s.kraken.report.txt" % (self.kraken_report_path, file.split(".")[0]), "w",
                  encoding="utf-8") as fw_kreOut:
            # STEP 1/4: READ TAXONOMY FILE
            count_nodes, root_node, taxId2node = 0, -1, {}

            with open("%s/ktaxonomy.tsv" % self.reference_database_path, 'r') as fr_taxFile:
                for line in fr_taxFile:
                    count_nodes += 1
                    [taxid, p_tid, rank, lvl_num, name] = line.strip().split('\t|\t')
                    curr_node = Tree(taxid, name, rank, lvl_num, p_tid)
                    taxId2node[taxid] = curr_node
                    # set parent/kids
                    if taxid == "1":
                        root_node = curr_node
                    else:
                        curr_node.parent = taxId2node[p_tid]
                        taxId2node[p_tid].add_child(curr_node)

            # STEP 2/4: READ KRAKEN FILE FOR COUNTS PER TAXID
            read_count, taxid2counts, taxid2allcounts = 0, {}, {}

            with open("%s/%s" % (self.kraken_output_path, file), 'r') as fr_kraken_file:
                for line in fr_kraken_file:
                    read_count += 1
                    l_vals = line.strip().split('\t')
                    taxid = l_vals[2]
                    # add to dictionaries
                    if taxid not in taxid2counts:
                        taxid2counts[taxid] = 1
                        taxid2allcounts[taxid] = 1
                    else:
                        taxid2counts[taxid] += 1
                        taxid2allcounts[taxid] += 1
            for curr_tid in taxid2counts:
                # Skip unclassified
                if curr_tid == '0':
                    continue
                p_node = taxId2node[curr_tid].parent
                add_counts = taxid2counts[curr_tid]
                # Assign reads for node
                taxId2node[curr_tid].lvl_reads += add_counts
                taxId2node[curr_tid].all_reads += add_counts
                while (p_node != None):
                    # Add child reads to parent node
                    p_taxid = p_node.tax_id
                    if p_taxid not in taxid2allcounts:
                        taxid2allcounts[p_taxid] = add_counts
                    else:
                        taxid2allcounts[p_taxid] += add_counts
                    p_node.all_reads += add_counts
                    # Get next parent node
                    p_node = p_node.parent
            # STEP 4/4: PRINT REPORT FILE
            # Write line for unclassified reads:
            if '0' in taxid2counts:
                fw_kreOut.write("%6.2f\t" % (float(taxid2counts['0']) / float(read_count) * 100))
                fw_kreOut.write("%i\t%i\t" % (taxid2counts['0'], taxid2counts['0']))
                fw_kreOut.write('U\t0\tunclassified\n')
            # Get remaining lines
            parse_nodes = [root_node]
            while len(parse_nodes) > 0:
                curr_node = parse_nodes.pop(0)
                curr_tid = curr_node.tax_id
                # Print information for this level
                fw_kreOut.write("%6.2f\t" % (float(taxid2allcounts[curr_tid]) / float(read_count) * 100))
                fw_kreOut.write("%i\t" % taxid2allcounts[curr_tid])
                if curr_tid not in taxid2counts:
                    fw_kreOut.write("0\t")
                else:
                    fw_kreOut.write("%i\t" % taxid2counts[curr_tid])
                fw_kreOut.write("%s\t" % curr_node.level_rank)
                fw_kreOut.write("%s\t" % curr_tid)
                fw_kreOut.write(" " * curr_node.level_num * 2 + curr_node.name + "\n")
                # Add children to list
                for child in sorted(curr_node.children, key=operator.attrgetter('all_reads')):
                    if child.tax_id not in taxid2allcounts:
                        continue
                    if taxid2allcounts[child.tax_id] == 0:
                        continue
                    # Add to list
                    parse_nodes.insert(0, child)

    @time_it
    @write_logger
    def kraken2_step4_single_kraken(self):
        """
        Make kraken report for intraAlignment_analysis cell
        :return:
        """
        # print(self.step_style % "Step4: make kraken report for intraAlignment_analysis cell")
        if self.test_mode:
            return
        with Pool(self.single_thread) as p:
            arg_list = os.listdir(self.kraken_output_path)
            list(tqdm(p.imap(self.make_kreport_report, arg_list), ncols=80, total=len(arg_list),
                      desc="Progressing"))

    # make bracken report for intraAlignment_analysis cell
    def make_bracken_report(self, file):
        """
        Make bracken report for intraAlignment_analysis cell
        :param file:
        :return:
        """
        command6 = "bracken -d %s -l S " \
                   "-i %s/%s.kraken.report.txt " \
                   "-w %s/%s.bracken.report.txt " \
                   "-o %s/%s.bracken.species.report.txt >/dev/null 2>&1" \
                   % (self.reference_database_path,
                      self.kraken_report_path, file.split(".")[0],
                      self.bracken_report_path, file.split(".")[0],
                      self.bracken_result_path, file.split(".")[0])
        subprocess.run(command6, shell=True)

    @time_it
    @write_logger
    def kraken2_step5_single_bracken(self):
        """
        Make bracken report for intraAlignment_analysis cell
        :return:
        """
        # print(self.step_style % "Step5: make bracken report for intraAlignment_analysis cell")
        if self.test_mode:
            return
        with Pool(self.single_thread) as p:
            arg_list = os.listdir(self.kraken_report_path)
            list(tqdm(p.imap(self.make_bracken_report, arg_list), ncols=80, total=len(arg_list),
                      desc="Progressing"))

    ########################################################
    # refer to "https://github.com/jenniferlu717/KrakenTools/blob/master/kreport2mpa.py"
    ########################################################
    def kraken2mpa(self, file):
        """
        Convert bracken report to mpa format
        :param file:
        :return:
        """

        def process_kraken_report(curr_str):
            """
            Process kraken report
            :param curr_str:
            :return:
            """
            split_str = curr_str.strip().split('\t')
            if len(split_str) < 4:
                return []
            try:
                int(split_str[1])
            except ValueError:
                return []
            percents = float(split_str[0])
            all_reads = int(split_str[1])
            # Extract relevant information
            try:
                taxid = int(split_str[-3])
                level_type = split_str[-2]
                map_kuniq = {'species': 'S', 'genus': 'G', 'family': 'F',
                             'order': 'O', 'class': 'C', 'phylum': 'P', 'superkingdom': 'D',
                             'kingdom': 'K'}
                if level_type not in map_kuniq:
                    level_type = '-'
                else:
                    level_type = map_kuniq[level_type]
            except ValueError:
                taxid = int(split_str[-2])
                level_type = split_str[-3]
            # Get name and spaces
            spaces = 0
            name = split_str[-1]
            for char in name:
                if char == ' ':
                    name = name[1:]
                    spaces += 1
                else:
                    break
            name = name.replace(' ', '_')
            # Determine level based on number of spaces
            level_num = spaces / 2
            return [name, level_num, level_type, all_reads, percents]

        # Process report file and output
        curr_path, prev_lvl_num = [], -1
        # Read through report file
        main_lvls = ['R', 'K', 'D', 'P', 'C', 'O', 'F', 'G', 'S']
        with open("%s/%s.bracken.report.txt" % (self.bracken_report_path, file.split(".")[0]), 'r') as fr_file:
            with open("%s/%s.bracken.mpa.txt" % (self.mpa_report_path, file.split(".")[0]), 'w') as fw_file:
                for line in fr_file:
                    report_vals = process_kraken_report(line)
                    # If header line, skip
                    if len(report_vals) < 5:
                        continue
                    # Get relevant information from the line
                    [name, level_num, level_type, all_reads, percents] = report_vals
                    if level_type == 'U':
                        continue
                    # Create level name
                    if level_type not in main_lvls:
                        level_type = "x"
                    elif level_type == "K":
                        level_type = "k"
                    elif level_type == "D":
                        level_type = "k"
                    level_str = level_type.lower() + "__" + name
                    # Determine full string to add
                    if prev_lvl_num == -1:
                        # First level
                        prev_lvl_num = level_num
                        curr_path.append(level_str)
                    else:
                        # Move back if needed
                        while level_num != (prev_lvl_num + 1):
                            prev_lvl_num -= 1
                            curr_path.pop()
                        # Print if at non-traditional level and that is requested
                        if level_type != "x":
                            # Print all ancestors of current level followed by |
                            for string in curr_path:
                                if string[0] != "x":
                                    if string[0] != "r":
                                        fw_file.write(string + "|")
                            # Print final level and then number of reads
                            fw_file.write(level_str + "\t" + str(all_reads) + "\n")
                        # Update
                        curr_path.append(level_str)
                        prev_lvl_num = level_num

    ########################################################
    # refer to "https://github.com/jenniferlu717/KrakenTools/blob/master/combine_mpa.py"
    ########################################################
    def combine_mpa(self):
        """
        Combine mpa report
        :return:
        """
        samples, sample_count, values, parent2child, toparse = {}, 0, {}, {}, []
        for in_file in os.listdir(self.mpa_report_path):
            with open("%s/%s" % (self.mpa_report_path, in_file), 'r') as fr_file:
                sample_count += 1
                sample_name = in_file.split(".")[0]
                for line in fr_file:
                    [classification, val] = line.strip().split('\t')
                    # Check for parents
                    split_vals = classification.split("|")
                    curr_parent = ''
                    for i in range(0, len(split_vals)):
                        test_val = "|".join(split_vals[0:i])
                        if test_val in values:
                            curr_parent = test_val
                            # No parent
                    if curr_parent == '':
                        if classification not in values:
                            toparse.append(classification)
                            # Most specific parent found
                    if curr_parent != '':
                        if curr_parent not in parent2child:
                            parent2child[curr_parent] = []
                        if classification not in parent2child[curr_parent]:
                            parent2child[curr_parent].append(classification)
                    # Save classification to value map
                    if classification not in values:
                        values[classification] = {}
                    values[classification][sample_count] = val
                # Save sample name
                samples[sample_count] = sample_name
        # Write header
        with gzip.open("%s/combinedMpa.kraken2.txt.gz" % self.output_path, 'wt') as fw_file:
            fw_file.write("#Classification")
            for i in range(1, sample_count + 1):
                fw_file.write("\t" + samples[i])
            fw_file.write("\n")

            # Write each line
            parsed = {}
            count_c = 0
            while len(toparse) > 0:
                curr_c = toparse.pop(0)
                # Add all children to stack
                if curr_c in parent2child:
                    for child in parent2child[curr_c]:
                        toparse.insert(0, child)
                        # For the current classification, print per sample
                fw_file.write(curr_c)
                for i in range(1, sample_count + 1):
                    if i in values[curr_c]:
                        fw_file.write("\t" + values[curr_c][i])
                    else:
                        fw_file.write("\t0")
                fw_file.write("\n")
                count_c += 1

    @time_it
    @write_logger
    def kraken2_step6_convert_to_map_format(self):
        """
        Convert bracken report to mpa format and combined output results
        :return:
        """
        # print(self.step_style % "Step6: convert bracken report to mpa format and combined output results")

        if self.test_mode:
            return
        with Pool(self.single_thread) as p:
            arg_list = os.listdir(self.bracken_report_path)
            list(tqdm(p.imap(self.kraken2mpa, arg_list), ncols=80, total=len(arg_list),
                      desc="Progressing"))
        self.combine_mpa()

        with gzip.open("%s/filtered_feature_bc_matrix/barcodes.tsv.gz" % self.ranger_out_path, "rt") as f:
            # check if the barcodes in the kraken2 output，if not, add the barcode to the output file
            barcodes = [x.strip() for x in f.readlines()]

            s_microbe = pd.read_csv("%s/combinedMpa.kraken2.txt.gz" % self.output_path, sep="\t", index_col=0,
                                    compression='gzip', encoding="utf-8")

            barcodes_nm = [x for x in barcodes if x not in s_microbe.columns]

            s_microbe = pd.concat([s_microbe,
                                   pd.DataFrame(data=np.zeros((s_microbe.shape[0], len(barcodes_nm))),
                                                columns=barcodes_nm,
                                                index=s_microbe.index)], axis=1)
            s_microbe = s_microbe.astype(int)
            s_microbe = s_microbe[barcodes]
            s_microbe.to_csv("%s/combinedMpa.kraken2.txt.gz" % self.output_path, sep="\t", compression='gzip',
                             encoding="utf-8")

    ########################################################
    # blast microbial alignment progress
    ########################################################
    @time_it
    @write_logger
    def blast_step1_extract_unmapped_reads(self):
        """
        Extract unmapped reads from the bam file
        :return:
        """

        # check test mode
        if self.test_mode:
            return

        # extract barcode information from bam file
        samtools_command = [
            "samtools view",
            "-@ 24",
            "-f 4",
            "%s/possorted_genome_bam.bam" % self.ranger_out_path
        ]
        awk_command = r'''awk 'match($0,/CB:Z:[ATCG]*-1/,b)&&match($0,/UB:Z:[ATCG]*/,u){print ">"$1"|"b[0]"|"u[0]"\n"$10}' ''' + \
                      '''> %s/possorted_genome_bam_unmapped_noHost.fasta''' % self.output_path

        # Combine the commands into a single string
        extract_command = " ".join(samtools_command) + " | " + awk_command
        subprocess.run(extract_command, shell=True, check=True)

    # blast alignment progress
    @time_it
    @write_logger
    def blast_step2_alignment(self):
        """
        blast alignment progress
        :return:
        """
        # check test mode
        if self.test_mode:
            return
        else:
            # blast alignment progress
            blast_command = ("blastn -query %s/possorted_genome_bam_unmapped_noHost.fasta "
                             "-db %s/nt "
                             "-qcov_hsp_perc 80 -perc_identity 80 "
                             "-outfmt '6 qseqid staxid stitle length evalue pident qcovus' "
                             "-max_target_seqs 1 -max_hsps 1 -culling_limit 1 "
                             "-num_threads %s "
                             "> %s/unmapped_noHost_blast.tsv") % (
                                self.output_path, self.reference_database_path,
                                self.blast_thread, self.output_path)
            subprocess.run(blast_command, shell=True)

    # convert blastn result to matrix
    @time_it
    @write_logger
    def blast_step3_convert_to_matrix(self):
        """
        Convert blastn result to matrix
        :return:
        """
        # split out the barcode and umi columns
        blast_command_1 = "cat %s/unmapped_noHost_blast.tsv | sed 's/|[CU]B:Z:/\t/g' | cut -f1-5 > %s/unmapped_noHost_blast.tsv.temp" % (
            self.output_path, self.output_path)
        subprocess.run(blast_command_1, shell=True)

        df = pd.read_table("%s/unmapped_noHost_blast.tsv.temp" % self.output_path, header=None, dtype=str)
        # Extract unique barcodes and taxids
        barcode_list = list(np.unique(df[1]))
        taxid_list = list(np.unique(df[3]))

        # Initialize the result DataFrame
        matrix = pd.DataFrame(index=taxid_list, columns=barcode_list, dtype='int64').fillna(0)

        def process_barcode(barcode_):
            """
            Process barcode
            :param barcode_:
            :return:
            """
            df_barcode = df[df[1] == barcode_]
            df_umi = [df_barcode[df_barcode[2] == x] for x in list(np.unique(df_barcode[2]))]

            barcode_taxid = [x[3].value_counts().index[0] if len(x[3].value_counts()) == 1 else
                             x[3].value_counts().index[0] if (x[3].value_counts().iloc[0] / len(x[3]) > 0.5 and
                                                              x[3].value_counts().iloc[0] != x[3].value_counts().iloc[
                                                                  1])
                             else "collision"
                             for x in df_umi]

            taxid_counts_ = pd.Series(barcode_taxid).value_counts()
            taxid_counts_ = taxid_counts_[taxid_counts_.index != "collision"]

            return barcode_, taxid_counts_

        # Use tqdm to add a progress bar
        for barcode in tqdm(barcode_list, ncols=80, desc="Processing barcodes"):
            barcode, taxid_counts = process_barcode(barcode)
            matrix.loc[taxid_counts.index, barcode] = taxid_counts.values

        # Save the result matrix to a file

        matrix.to_csv("%s/unmapped_noHost_blast_matrix.tsv.gz" % self.output_path, compression='gzip', encoding="utf-8",
                      sep="\t")

    @time_it
    @write_logger
    def blast_step4_convert_taxid_name(self):
        """
        Convert taxid to taxon name
        :return:
        """
        # paste <(echo "" && tail -n +2 matrix.mdx | cut -f 1 | taxonkit lineage --data-dir $DB| taxonkit reformat --data-dir $DB -F -f "k__{k}|p__{p}|c__{c}|o__{o}|f__{f}|g__{g}|s__{s}" | cut -f 3) <(cut --complement -f 1 matrix.mdx) > matrix.mpa.mdx.temp
        subprocess.run(
            ["zcat %s/unmapped_noHost_blast_matrix.tsv.gz > %s/unmapped_noHost_blast_matrix.tsv" % (self.output_path,
                                                                                                    self.output_path)],
            shell=True, check=True)
        blast_command_2 = ("paste <(echo '' && tail -n +2 %s/unmapped_noHost_blast_matrix.tsv "
                           "| cut -f 1 | taxonkit lineage --data-dir %s "
                           "| taxonkit reformat --data-dir %s -F -f 'k__{k}|p__{p}|c__{c}|o__{o}|f__{f}|g__{g}|s__{s}' "
                           "| cut -f 3) <(cut --complement -f 1 %s/unmapped_noHost_blast_matrix.tsv) "
                           "> %s/unmapped_noHost_blast_matrix.mpa.tsv") % (
                              self.output_path, self.reference_database_path, self.reference_database_path,
                              self.output_path, self.output_path)
        print(blast_command_2)
        subprocess.run(["bash", "-c", blast_command_2], check=True)

    @time_it
    @write_logger
    def blast_step5_convert_to_map_format(self):
        """
        Convert blastn result to mpa format
        :return:
        """
        # Add data with the same row name
        with gzip.open("%s/filtered_feature_bc_matrix/barcodes.tsv.gz" % self.ranger_out_path, "rt") as f:
            # check if the barcodes in the kraken2 output，if not, add the barcode to the output file
            barcodes = [x.strip() for x in f.readlines()]

            s_microbe = pd.read_csv("%s/unmapped_noHost_blast_matrix.mpa.tsv" % self.output_path, sep="\t",
                                    index_col=0,
                                    encoding="utf-8")
            s_microbe = s_microbe.groupby(s_microbe.index).sum()
            barcodes_nm = [x for x in barcodes if x not in s_microbe.columns]

            s_microbe = pd.concat([s_microbe,
                                   pd.DataFrame(data=np.zeros((s_microbe.shape[0], len(barcodes_nm))),
                                                columns=barcodes_nm,
                                                index=s_microbe.index)], axis=1)
            s_microbe = s_microbe.astype(int)
            s_microbe = s_microbe[barcodes]
            # remove the row with all zeros
            s_microbe = s_microbe.loc[s_microbe.sum(axis=1) != 0]
            s_microbe.to_csv("%s/unmapped_noHost_blast_matrix.mpa.tsv" % self.output_path,
                             sep="\t",
                             encoding="utf-8")

        # { head -n 1 matrix.mpa.mdx.temp; tail -n +2 matrix.mpa.mdx.temp | grep '^k__Bacteria'; } > matrix.mpa.mdx
        blast_command_3 = ("{ head -n 1 %s/unmapped_noHost_blast_matrix.mpa.tsv; "
                           "tail -n +2 %s/unmapped_noHost_blast_matrix.mpa.tsv | grep '^k__Bacteria'; } "
                           "> %s/combinedMpa.blast.txt") % (
                              self.output_path, self.output_path, self.output_path)
        subprocess.run(blast_command_3, shell=True)

        blast_command_4 = "gzip %s/combinedMpa.blast.txt" % self.output_path
        subprocess.run(blast_command_4, shell=True)

    @time_it
    def run(self):
        """
        Run the pipeline
        :return:
        """
        if self.align_method == "kraken2":
            if self.platform == "10x":
                self.kraken2_step1_extract_unmapped_reads()
            else:
                self.c4_extract_unmapped_reads()
            self.kraken2_step2_alignment()
            self.kraken2_step3_split_kraken_output()
            self.kraken2_step4_single_kraken()
            self.kraken2_step5_single_bracken()
            self.kraken2_step6_convert_to_map_format()

            if self.clean_mode:
                # remove temporary files, keep the final results
                # directory: singleCell_kraken_output, singleCell_kraken_report, singleCell_bracken_report,
                # singleCell_bracken_result, singleCell_mpa_report
                paths_to_remove = [
                    "singleCell_kraken_output",
                    "singleCell_kraken_report",
                    "singleCell_bracken_report",
                    "singleCell_bracken_result",
                    "singleCell_mpa_report",
                    "unmapped_noHost_kraken2_filtered.output",
                    "kraken.report.txt"
                ]
                remove_command = "rm -rf " + " ".join(os.path.join(self.output_path, path) for path in paths_to_remove)

                subprocess.run(remove_command, shell=True)

        elif self.align_method == "blast":
            if self.platform == "10x":
                self.blast_step1_extract_unmapped_reads()
            else:
                self.c4_extract_unmapped_reads()
            self.blast_step2_alignment()
            self.blast_step3_convert_to_matrix()
            self.blast_step4_convert_taxid_name()
            self.blast_step5_convert_to_map_format()

            if self.clean_mode:
                # remove temporary files, keep the final results
                # files: unmapped_noHost_blast.tsv, unmapped_noHost_blast_matrix.tsv, unmapped_noHost_blast_matrix.mpa.tsv
                paths_to_remove = [
                    "unmapped_noHost_blast.tsv.temp",
                    "unmapped_noHost_blast_matrix.tsv",
                    "unmapped_noHost_blast_matrix.mpa.tsv"
                ]
                remove_command = "rm -rf " + " ".join(os.path.join(self.output_path, path) for path in paths_to_remove)
                subprocess.run(remove_command, shell=True)


def main():
    """
    Main function
    :return:
    """
    # Command line parameter
    parser = ap.ArgumentParser(
        description='Single cell microbiome alignment tool'
    )

    # Required arguments
    parser.add_argument(
        '-rp',
        '--rgPath',
        dest='ranger_out_path',
        help="Path of the 'outputs' directory generated by cellranger or spaceranger pipeline"
    )

    #
    parser.add_argument(
        "-method",
        "--method",
        dest="align_method",
        help="Alignment method, default: kraken2, choose from kraken2 and blast",
        choices=['kraken2', 'blast'],
        default="kraken2"
    )

    # platform choosing from 10 and c4, c4: DNBelab C4
    parser.add_argument(
        '-p',
        '--platform',
        dest='platform',
        default="10x",
        type=str,
        choices=["10x", "c4"],
        help='platform choosing from 10 and c4, 10x: 10x Genomics; c4: DNBelab C4, default: 10x'
    )

    parser.add_argument(
        '-c4_o', '--c4_output_path',
        dest='c4_output_path',
        default=None,
        type=str,
        help='C4 output path. If platform is c4, this path is required and must contain "01.data" and "output" folders.'
    )

    parser.add_argument(
        '-c4_r1', '--c4_cdna_r1',
        dest='c4_cdna_r1',
        default=None,
        type=str,
        help='Absolute path of R1 sequencing file of cDNA library. '
             'Required if platform is c4.  '
             'Filename must contain "cDNA_R1.fq"'
    )

    default_options = {
        'kraken_thread': 24,
        'single_thread': 20,
        'blast_thread': 24
    }

    parser.add_argument(
        '-db',
        '--dbPath',
        dest='reference_database_path',
        help="Kraken database path, if align method is 'kraken2', kraken2 reference database path is required, "
             "recommend download from 'https://benlangmead.github.io/aws-indexes/k2'; "
             "else if align method is 'blast',  blast reference database path is required",
        required=True
    )

    parser.add_argument(
        '-kt',
        '--kThread',
        dest='kraken2_thread',
        default=default_options["kraken_thread"],
        type=int,
        help='kraken2 process thread number, default: %s' % default_options["kraken_thread"]
    )

    parser.add_argument(
        '-st',
        '--sThread',
        dest='single_thread',
        default=default_options["single_thread"],
        type=int,
        help='intraAlignment_analysis cell process thread number, default: %s' % default_options["single_thread"]
    )

    parser.add_argument(
        '-bt',
        '--bThread',
        dest='blast_thread',
        default=default_options['blast_thread'],
        type=int,
        help='blast process thread number, default: %s' % default_options['blast_thread']
    )

    # Test mode, test mode will skip the alignment process, test the pipeline only, default: False.
    parser.add_argument(
        '-test',
        '--test_mode',
        dest='test_mode',
        action='store_true',
        default=False,
        help='test mode, if set to True, skip the alignment process, test the pipeline only, default: False'
    )

    # Clean temporary files during the process, default: False.
    parser.add_argument(
        '-clean',
        '--clean_mode',
        dest='clean_mode',
        action='store_true',
        default=False,
        help='clean mode, if set to True, clean temporary files during the process, default: False'
    )

    args = parser.parse_args()
    if args.platform == 'c4':
        if not args.c4_output_path:
            print("Error: --c4_output_path is required when platform is 'c4'.")
            sys.exit(1)

        if not os.path.isdir(args.c4_output_path):
            print(f"Error: {args.c4_output_path} is not a valid directory.")
            sys.exit(1)

        expected_subdirs = ['01.data', 'output']
        for subdir in expected_subdirs:
            subdir_path = os.path.join(args.c4_output_path, subdir)
            if not os.path.isdir(subdir_path):
                print(f"Error: {subdir_path} is missing. The c4_output_path must contain '{subdir}' folder.")
                sys.exit(1)

        # check c4_cdna_r1
        if not args.c4_cdna_r1:
            print("Error: '-c4_r1' or '--c4_cdna_r1' is required when platform is 'c4'.")
            sys.exit(1)
        if not os.path.isfile(args.c4_cdna_r1):
            print(f"Error: {args.c4_cdna_r1} is not a valid file.")
            sys.exit(1)
        if "cDNA_R1" not in os.path.basename(args.c4_cdna_r1):
            print(f"Error: {args.c4_cdna_r1} should be R1 sequencing file of cDNA library, "
                  f"and the filename should contain 'cDNA_R1'.")
            sys.exit(1)

    task = intraAlignment(ranger_out_path=args.ranger_out_path,
                          align_method=args.align_method,
                          reference_database_path=args.reference_database_path,
                          platform=args.platform,
                          c4_output_path=args.c4_output_path,
                          c4_cdna_r1=args.c4_cdna_r1,
                          kraken_thread=args.kraken2_thread,
                          single_thread=args.single_thread,
                          blast_thread=args.blast_thread,
                          test_mode=args.test_mode,
                          clean_mode=args.clean_mode)
    task.run()


if __name__ == '__main__':
    try:
        print(Fore.LIGHTBLUE_EX
              + "%-10s : %s" % ("start time", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        main()
        print(Fore.LIGHTBLUE_EX
              + "%-10s : %s" % ("end time", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    except Exception as e:
        print(Fore.LIGHTRED_EX + "Process terminated with errors: \n%s" % e)
        print(Fore.LIGHTBLUE_EX
              + "%-10s : %s" % ("end time", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

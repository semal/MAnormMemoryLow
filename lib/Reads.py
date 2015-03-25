"""
this module is used for processing reads, including reads reading and splitting.
"""
import os


class Reads:
    """
    Reads class defines the data type of rds reading from rds file
    """

    def __init__(self):
        """
        rds including marker_name, tag_shift, tag_len, chr_id, position, valuence, six properties.
        """
        self.mk_name = ''
        self.chr_id = []
        self.pos = []   # reads position

    @staticmethod
    def read_rds_file(rds_file_path, tag_shift, tag_len):
        """
        read reads data from a reads file, including non_split reads file and split reads files
        @param rds_file_path: reads file path
        @param tag_shift:
            shift size of reads, which equal to the half of the DNA fragment after size selection.
        @param tag_len:
            length of reads
        """
        rds = Reads()
        with open(rds_file_path) as f:
            for line in f:
                split_line = line.strip().split()
                if split_line[-1] == '+':
                    rds.pos.append(int(split_line[1]) + tag_shift)
                elif split_line[-1] == '-':
                    rds.pos.append(int(split_line[2]) - tag_shift)
        return rds

    @staticmethod
    def split_rds_by_chr(rds_file_path, chr_id_list, output_dir):
        """
        Because rds file generally is much big, so splitting it by each chromosome into each file
        @param rds_file_path: reads file path
        @param chr_id_list: chromosomes id list for splitting reads by chromosome
        @param output_dir: folder of split reads saving
        """
        # create folder of split rds data files
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        folder_of_split_data = os.sep.join([output_dir, 'split_rds_data'])
        if not os.path.exists(folder_of_split_data):
            os.mkdir(folder_of_split_data)

        # create subfolder and files of each chromosome
        all_f = {}
        marker_name = os.path.basename(rds_file_path)[:-4]
        for chr_id in chr_id_list:
            subfolder_of_chr = os.sep.join([folder_of_split_data, chr_id])
            if not os.path.exists(subfolder_of_chr):
                os.mkdir(subfolder_of_chr)
            f = open(os.sep.join([subfolder_of_chr, '_'.join([marker_name, chr_id + '.txt'])]), 'w')
            all_f.update({chr_id.lower(): f})

        # split data into each file
        rds_f = open(rds_file_path)
        for line in rds_f:
            split_line = line.split()
            if split_line[0].lower() in all_f.keys():
                all_f[split_line[0].lower()].write(line)

        # close all split files
        for chr_id in chr_id_list:
            all_f[chr_id.lower()].close()

    @staticmethod
    def read_split_rds_file(marker_name, chr_id, tag_shift, tag_len, output_dir):
        """
        read split rds through chr_id
        @param marker_name: name of the reads file
        @param chr_id: chromosome name
        @param tag_shift: tag shift size
        @param tag_len: read length
        @param output_dir: the folder of reads split
        """
        rds_file_path = os.sep.join([output_dir, 'split_rds_data', chr_id, '_'.join([marker_name, chr_id + '.txt'])])
        rds = Reads.read_rds_file(rds_file_path, tag_shift, tag_len)
        rds.mk_name = marker_name
        return rds
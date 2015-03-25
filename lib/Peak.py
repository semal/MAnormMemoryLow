"""
this module is used for processing peaks, including peaks reading, choosing, sorting, and reads density calculating.
"""
from bisect import bisect_left, bisect_right
import os

import numpy as np
from numpy import array

import Reads


class Peak:
    """
    Define a peak class to conveniently get and store macs xls peak's traits
    """
    def __init__(self):
        """
        peak including name, pk_num, chr_id, start, end, summit, six properties.
        But another one added after calculating the peak reads density.
        """
        self.name = ''
        self.pk_num = 0
        self.chr_id = []
        self.start = []
        self.end = []
        self.summit = []
        self.rds_count = {}
        self.rds_density = {}
        self.another_info = {}

    @staticmethod
    def read_pk_file(pk_file_path):
        """
        read peak from a peak file
        @param pk_file_path: peak file path
        """
        pk = Peak()
        pk.name = os.path.basename(pk_file_path)[:-4]
        if pk.name.endswith('_peaks'):
            pk.name = pk.name.replace('_peaks', '')

        pk_file = open(pk_file_path)
        summit_column = 0
        for line in pk_file.readlines():
            # get the column_number of summit
            if line.lower().startswith('chr') and len(line.split()[0]) == 3:
                split_line = array(line.split())
                idx_of_summit_column = np.where(split_line == 'summit')[0]
                if idx_of_summit_column.size == 1:
                    summit_column = idx_of_summit_column[0]

            if line.lower().startswith('chr') and len(line.split()[0]) > 3:
                split_line = line.split()
                pk.chr_id.append(split_line[0])
                pk.start.append(int(split_line[1]))
                pk.end.append(int(split_line[2]))
                if not summit_column == 0:
                    pk.summit.append(int(split_line[summit_column]) + int(split_line[1]))
                else:
                    summit = (int(split_line[1]) + int(split_line[2])) / 2.0
                    pk.summit.append(int(round(summit)))
        pk_file.close()
        pk.pk_num = len(pk.chr_id)
        return pk

    def add(self, pk):
        """
        peak add peak
        @param pk: another Peak instance
        """
        addition_pk = Peak()
        addition_pk.name = self.name
        addition_pk.pk_num = self.pk_num + pk.pk_num
        addition_pk.chr_id = self.chr_id + pk.chr_id
        addition_pk.start = self.start + pk.start
        addition_pk.end = self.end + pk.end
        addition_pk.summit = self.summit + pk.summit
        for key in self.rds_density.keys():
            addition_pk.rds_density.update({key: array(self.rds_density[key].tolist() + pk.rds_density[key].tolist())})
        for key in self.another_info.keys():
            addition_pk.another_info.update(
                {key: array(self.another_info[key].tolist() + pk.another_info[key].tolist())})

        return addition_pk

    def choose_pks(self, idx_array):
        """
        choose peaks from a known peak with a array of idxes
        @param idx_array: index array for choosing peaks by this order
        """
        chose_pk = Peak()
        chose_pk.name = self.name
        chose_pk.pk_num = idx_array.size
        chose_pk.chr_id = array(self.chr_id)[idx_array].tolist()
        chose_pk.start = array(self.start)[idx_array].tolist()
        chose_pk.end = array(self.end)[idx_array].tolist()
        chose_pk.summit = array(self.summit)[idx_array].tolist()
        for key in self.rds_density.keys():
            chose_pk.rds_density.update({key: self.rds_density[key][idx_array]})
        for key in self.another_info.keys():
            chose_pk.another_info.update({key: self.another_info[key][idx_array]})

        return chose_pk

    def choose_one_chr_pks(self, chr_id):
        """
        choose peaks from a known peak with the same chr_id
        @param chr_id: chromosome name
        """
        idxes = np.where(array(self.chr_id) == chr_id)[0]
        pk = self.choose_pks(idxes)
        return pk, idxes

    def sort_pks(self, value):
        """
        sort peaks by summits or starts
        @param: peak summit or start point used for sorting peaks.
        """
        pk_sorted = Peak()
        if value == 'summit':
            sorted_idxes = np.argsort(array(self.summit))
        elif value == 'start':
            sorted_idxes = np.argsort(array(self.start))
        else:
            raise ValueError()
        pk_sorted.name = self.name
        pk_sorted.pk_num = self.pk_num
        for idx in sorted_idxes.tolist():
            pk_sorted.chr_id.append(self.chr_id)
            pk_sorted.start.append(self.start[idx])
            pk_sorted.end.append(self.end[idx])
            pk_sorted.summit.append(self.summit[idx])
        for key in self.rds_density.keys():
            density = []
            for idx in sorted_idxes.tolist():
                density.append(self.rds_density[key][idx])
            pk_sorted.rds_density.update({key: array(density)})
        for key in self.another_info.keys():
            val = []
            for idx in sorted_idxes.tolist():
                val.append(self.another_info[key][idx])
            pk_sorted.another_info.update({key: array(val)})

        return pk_sorted

    def cal_rds_density(self, mk_name, ext, tag_shift, tag_len, output_dir):
        """
        using reads data to calculate peak read density of this marker
        @param mk_name: reads name
        @param ext: extension size, used to extent peak length
        @param tag_shift: tag shift
        @param tag_len: tag length
        @param output_dir: folder of split reads files saving
        """
        if mk_name not in self.rds_density.keys():
            print 'calculating %s read density using %s reads' % (self.name, mk_name)
            rds_density = np.zeros(self.pk_num)
            rds_count = np.zeros(self.pk_num)
            chr_id_list = list(set(self.chr_id))
            for chr_id in chr_id_list:
                print '%s ...' % chr_id
                pk_chr, idxes_chr = self.choose_one_chr_pks(chr_id)
                rds = Reads.read_split_rds_file(mk_name, chr_id, tag_shift, tag_len, output_dir)
                rds_pos = array(rds.pos)
                rds_pos.sort()

                density_chr = []
                count_chr = []
                for idx in xrange(pk_chr.pk_num):
                    start = pk_chr.summit[idx] - ext - 1
                    end = pk_chr.summit[idx] + ext

                    # rd_count = rds_pos[(rds_pos <= end) & (rds_pos >= start)].size
                    rd_count = cal_read_count(start, end, rds_pos)
                    one_density = (rd_count + 1) * 1000 / (2.0 * ext)

                    density_chr.append(one_density)
                    count_chr.append(rd_count)
                rds_density[idxes_chr] = array(density_chr)
                rds_count[idxes_chr] = array(count_chr)

            self.rds_density.update({mk_name: rds_density})
            self.rds_count.update({mk_name: rds_count})


def cal_read_count(peak_start, peak_end, reads_position):
    # this idea is from Limushan
    si = bisect_left(reads_position, peak_start)
    ei = bisect_right(reads_position, peak_end)
    try:
        if peak_end == reads_position[ei]:
            return ei - si + 1
        else:
            return ei - si
    except IndexError:
        return ei - si


if __name__ == '__main__':
    pass
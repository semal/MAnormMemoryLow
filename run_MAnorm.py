import os
import time


start = time.time()
command = 'python MAnorm ' \
          '--p1 1_peaks ' \
          '--r1 1-reads.bed ' \
          '--l1 36 ' \
          '--s1 100 ' \
          '--p2 2_peaks ' \
          '--r2 2-reads.bed ' \
          '--l2 36 ' \
          '--s2 100 ' \
          '-n 2 ' \
          '-e 1000 ' \
          '-s ' \
          '--output test'
print command
# os.system(command)
print time.time() - start
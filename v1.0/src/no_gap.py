import gzip
import shutil
import re
from dataclasses import dataclass, field, asdict
from typing import List
import os
import pandas as pd
from concurrent import futures
import time
import datetime
# Download the transcript
#!wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.lncRNA_transcripts.fa.gz

with gzip.open('./gencode.v44.pc_transcripts.fa.gz', 'rb') as file_in:
    with open('./gencode.v44.pc_transcripts.fa', 'wb') as file_out:
        shutil.copyfileobj(file_in, file_out)
        print('fasta file created')

@dataclass
class Sequence:
    seq: str
    start: int
    end: int

@dataclass
class Stem:
    loop: Sequence
    upstream: Sequence = None
    downstream: Sequence = None

@dataclass
class Triplex:
    H: Sequence
    C: Sequence  # Swapped from W to C
    W: Sequence  # Swapped from C to W
    upperStem: Stem
    lowerStem: Stem
    pseudoknot_seq: Sequence = None
    pseudoknot_status: str = "NO"

    def to_dict(self):
        d = {
            'Hoogsteen Strand': self.H.seq,
            'H Start Index': self.H.start + 1,
            'H End Index': self.H.end + 1,

            'Crick Strand': self.C.seq,
            'C Start Index': self.C.start + 1,
            'C End Index': self.C.end + 1,

            'Watson Strand': self.W.seq,
            'W Start Index': self.W.start + 1,
            'W End Index': self.W.end + 1,

            'Upper Loop Sequence': self.upperStem.loop.seq,
            'Upper Loop Start Index': self.upperStem.loop.start + 1,
            'Upper Loop End Index': self.upperStem.loop.end + 1,

            'Lower Loop Sequence': self.lowerStem.loop.seq,
            'Lower Loop Start Index': self.lowerStem.loop.start + 1,
            'Lower Loop End Index': self.lowerStem.loop.end + 1,

            'Lower Stem (downstream) Sequence': None if not self.lowerStem.downstream else self.lowerStem.downstream.seq,
            'Lower Stem (downstream) Start Index': None if not self.lowerStem.downstream else self.lowerStem.downstream.start + 1,
            'Lower Stem (downstream) End Index': None if not self.lowerStem.downstream else self.lowerStem.downstream.end + 1,

            'Lower Stem (upstream) Sequence': None if not self.lowerStem.upstream else self.lowerStem.upstream.seq,
            'Lower Stem (upstream) Start Index': None if not self.lowerStem.upstream else self.lowerStem.upstream.start + 1,
            'Lower Stem (upstream) End Index': None if not self.lowerStem.upstream else self.lowerStem.upstream.end + 1,

            # No bulge for now----

            'Upper Stem (downstream) Sequence': None if not self.upperStem.downstream else self.upperStem.downstream.seq,
            'Upper Stem (downstream) Start Index': None if not self.upperStem.downstream else self.upperStem.downstream.start + 1,
            'Upper Stem (downstream) End Index': None if not self.upperStem.downstream else self.upperStem.downstream.end + 1,

            'Upper Stem (upstream) Sequence': None if not self.upperStem.upstream else self.upperStem.upstream.seq,
            'Upper Stem (upstream) Start Index': None if not self.upperStem.upstream else self.upperStem.upstream.start + 1,
            'Upper Stem (upstream) End Index': None if not self.upperStem.upstream else self.upperStem.upstream.end + 1,

            'Pseudoknot Status': self.pseudoknot_status,
            'Pseudoknot Sequence': None if not self.pseudoknot_seq else self.pseudoknot_seq.seq,
            'Pseudoknot Start Index': None if not self.pseudoknot_seq else self.pseudoknot_seq.start + 1,
            'Pseudoknot End Index': None if not self.pseudoknot_seq else self.pseudoknot_seq.end + 1,
        }

        return d


from typing import List

covalent_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


def generate_patterns(s: str) -> list:
  """generates sequence with single gap"""
  res = []
  for i in range(3, len(s)-3):
    res.append(s[:i] + '.' + s[i:])
  return res


def build_comp_seq(s: str):
  return ''.join(covalent_map[c] for c in s[::-1])


def find_upper_stem(seq: str, start: int, end: int):

  for size in range((end-start)//2, 2, -1):
    upstream = seq[start: start+size]
    downstream = build_comp_seq(upstream)

    if downstream == seq[end-size: end]:
      return [(start, start+size), (end-size, end)]



def find_lower_stem(seq: str, hstart: int, wend: int, cstart):

  for size in range((cstart-wend)//2, 2, -1):
    upstream = seq[hstart-size: hstart]
    downstream = build_comp_seq(upstream)

    try:
      comp_idx = seq[wend: cstart].index(downstream)
      comp_idx += wend
    except Exception:
      comp_idx = -1

    if comp_idx > 0:
      return [(hstart-size, hstart), (comp_idx, comp_idx+size)]



def detect_triplex_seq(seq: str) -> List[Triplex]:
    seq_length = len(seq)
    seq = seq.upper().replace('T', 'U')

    res = {}
    start = 4
    #max_bulge_length = 10  # Adjust as needed
    extended_c_search_length = 250  # Adjust as needed

    while start < seq_length-200:
        for size in range(15, 6, -1):
            end = start+size
            H = seq[start: end]

            patterns = []

            #patterns = generate_patterns(H[::-1])
            patterns.append(H[::-1])
            for pattern in patterns:
                for match in re.finditer(pattern, seq[end: end+len(H)+101]):
                    W = match.group(0)
                    C = build_comp_seq(W)

                    cend = end + match.end()
                    cmatch = re.search(C, seq[cend: cend+len(C)+100])
                    if cmatch:
                        hseq = Sequence(pattern[::-1], start, end)
                        wseq = Sequence(W, end+match.start(), end+match.end())
                        cseq = Sequence(C, cend+cmatch.start(), cend+cmatch.end())
                        upper_stem = lower_stem = None

                        comp_seq = find_upper_stem(seq, hseq.end, wseq.start)
                        if comp_seq:
                            uidx, didx = comp_seq
                            upstream = Sequence(seq[uidx[0]: uidx[1]], uidx[0], uidx[1])
                            downstream = Sequence(seq[didx[0]: didx[1]], didx[0], didx[1])
                            loop = Sequence(seq[uidx[1]: didx[0]], uidx[1], didx[0])
                            if abs(uidx[1] - didx[0]) >= 3:
                                upper_stem = Stem(loop, upstream, downstream)
                                print(f"Upstream: {upstream.seq}, Downstream: {downstream.seq}") # Debugging line
                        else:
                            loop = Sequence(seq[hseq.end: wseq.start], hseq.end, wseq.start)
                            if abs(wseq.start - hseq.end) >= 3:
                                upper_stem = Stem(loop)

                        comp_seq = find_lower_stem(seq, hseq.start, wseq.end, cseq.start)
                        if comp_seq:
                            uidx, didx = comp_seq
                            upstream = Sequence(seq[uidx[0]: uidx[1]], uidx[0], uidx[1])
                            downstream = Sequence(seq[didx[0]: didx[1]], didx[0], didx[1])
                            loop = Sequence(seq[didx[1]: cseq.start], didx[1], cseq.start)
                            if abs(didx[1] - cseq.start) >= 3:
                                lower_stem = Stem(loop, upstream, downstream)
                        else:
                            loop = Sequence(seq[wseq.end: cseq.start], wseq.end, cseq.start)
                            if abs(wseq.end - cseq.start) >= 3:
                                lower_stem = Stem(loop)

                        if upper_stem:
                            upper_loop_seq = upper_stem.loop.seq
                        else:
                            upper_loop_seq = seq[hseq.end: wseq.start]

                        comp_upper_loop_seq = build_comp_seq(upper_loop_seq)
                        extended_c_end = cend + cmatch.end() + extended_c_search_length
                        extended_C = seq[cend: extended_c_end]

                        pseudoknot_start = extended_C.find(comp_upper_loop_seq)
                        if pseudoknot_start != -1:
                            pseudoknot_end = pseudoknot_start + len(comp_upper_loop_seq)
                            pseudoknot_seq = Sequence(extended_C[pseudoknot_start:pseudoknot_end], cend + pseudoknot_start, cend + pseudoknot_end)
                            #bulge_seq = None
                            #if pseudoknot_start > 0 and pseudoknot_start <= max_bulge_length:
                                #bulge_seq = Sequence(extended_C[:pseudoknot_start], cend, cend + pseudoknot_start)
                            if upper_stem and lower_stem:
                              res[H] = Triplex(hseq, wseq, cseq, upper_stem, lower_stem, pseudoknot_seq, "YES")
                        else:
                            if upper_stem and lower_stem:
                              res[H] = Triplex(hseq, wseq, cseq, upper_stem, lower_stem, None, "NO")

            if H in res:
                start = end
                break
        start += 1

    return list(res.values())


# with open('./NONHSAT229874.1.fa') as f:
#   f.readline()
#   seq = f.readline()
#
# # print(seq)
# # print(len(seq))
# print(seq[5524: 5529].replace('T', 'U'))
# print(seq[5584: 5589].replace('T', 'U'))
#
# res = detect_triplex_seq(seq)
#
# for tpx in res:
#   print(asdict(tpx))
#   print('H - ', tpx.H)
#   print('C - ', tpx.C)
#   print('W - ', tpx.W)
#   print()
#   print('Upper Stem(upstream) - ', tpx.upperStem.upstream)
#   print('Upper Stem(downstream) - ', tpx.upperStem.downstream)
#   print('Upper Stem(loop) - ', tpx.upperStem.loop)
#   print()
#   print('Lower Stem(upstream) - ', tpx.lowerStem.upstream)
#   print('Lower Stem(downstream) - ', tpx.lowerStem.downstream)
#   print('Lower Stem(loop) - ', tpx.lowerStem.loop)
# # generate_patterns('AAUUUGGUC')

"""**Reading fast files for all transcripts**"""



import os
import pandas as pd


working_dir = './'
file_items = sorted(os.listdir(working_dir))
gene_items = []


for file in file_items:
  if file.endswith('.fa'):
    try:
      with open(os.path.join(working_dir, file)) as f:
        file = None
        seq = ''

        while True:
          line = f.readline().replace('\n', '')
          if line.startswith('>'):
            # marks the end of previous sequence
            if file and seq:
              #if file.split('|')[-3] =='XACT':
              gene_items.append((file, seq))

            # starting another sequence
            file, seq = line, ''
          else:
            seq += line

          if not line:
            break

        if file and seq:
          #if file.split('|')[-3] =='XACT':
          gene_items.append((file, seq))

    except:
      print('Exception while reading fasta sequence in file.')

print(f'Collected {len(gene_items)} gene items') # To check the length of the total gene items
print(gene_items[0]) # To check the nature of the first collected gene items


"""## Running in multiprocessing"""

import time
import datetime
from concurrent import futures

CORES = 18


def process_seq(item):
  file, seq = item
  print(f'Processing {file}')
  res = detect_triplex_seq(seq)
  if res:
    print(f'\tFound {len(res)} result for file: {file}')

  rows = []
  for tpx in res:
    # print(tpx)
    row = {'filename': file}
    row.update(tpx.to_dict())
    rows.append(row)
  return rows


if __name__ == '__main__':
  print(f'Started At: {datetime.datetime.now()}')
  s = time.perf_counter()
  print(f'Starting for {len(gene_items)} gene items')

  # starting multi-processing...
  csv_rows = []
  with futures.ProcessPoolExecutor(max_workers=CORES) as exec:
    results = exec.map(process_seq, gene_items[:]) # Processing all gene items from all collected gene items
    for result in results:
      csv_rows.extend(result)


  # post multi-processing...
  e = time.perf_counter()
  print(f'Took {int((e-s)//60)} mins {int((e-s)%60)} seconds')

  end_date = datetime.datetime.now()
  end_date = end_date.replace(microsecond=0).isoformat()
  print(f'Ended At: {end_date}')

  # writing to csv
  df = pd.DataFrame(csv_rows)
  df.to_csv(f'triplex_pseudoknot_no_pseudo_buldge{end_date}_all_no_gap_minsize_7_pcRNA.csv', index=False)
  print(df)

  # for i, (f, s) in enumerate(gene_items):
  #   if f == '>ENST00000667414.1|ENSG00000228794.12|OTTHUMG00000196004.1|OTTHUMT00000523925.1|LINC01128-223|LINC01128|1752|':
  #     print(i)
  #
  # process_seq(gene_items[94])

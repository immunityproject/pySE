import csv
from collections import defaultdict

lanl = defaultdict(dict)
with open('lanl_ctl_summary.csv', 'rb') as lanlfile:
  reader = csv.reader(lanlfile)
  for row in reader:
    if len(row) < 8:
      print('Short row: {}'.format(row))
      continue
    if row[0] == 'Epitope':
      continue
    key = '{},{},{},{}'.format(row[0], row[1].upper(), row[2], row[3])
    if key in lanl:
      print('Duplicate key: {}'.format(key))
      continue
    lanl[key] = { 'peptide': row[0], 'protein': row[1].upper(), 'start': row[2], 'end': row[3], 'subprotein': row[4].upper(), 'subtype': row[6].upper() }

ipeps = []
with open('immunityproject-hiv-epitopes.csv', 'rb') as ipfile:
  reader = csv.reader(ipfile)
  for row in reader: 
    if len(row) < 5:
      continue
    if row[0] == 'Sequence':
      continue
    ipeps.append({ 'peptide': row[0], 'protein': row[1].upper(), 'range': row[2], 'epitope': row[3].upper() })

eps = []
for e in ipeps:
  start, end = e['range'].split('-')
  protein = e['protein']
  if protein == 'ENV':
    protein = 'GP160'
  lkey = '{},{},{},{}'.format(e['peptide'], protein, start, end)
  lanl_e = lanl[lkey]
  subprotein = lanl_e.get('subprotein', '')

  new_e = { k: v for k, v in e.items() if k not in ['range'] }
  if subprotein:
    subprotein = subprotein.split('(')[0]
  new_e['subprotein'] = subprotein
  new_e['start'] = start
  new_e['end'] = end
  new_e['subtype'] = lanl_e.get('subtype', '')
  eps.append(new_e)

with open('epitope-info.csv', 'wb') as out:
  fns = ['peptide', 'protein', 'subprotein', 'start', 'end', 'epitope', 'subtype']
  writer = csv.DictWriter(out, fieldnames=fns)
  writer.writeheader()
  for e in eps:
    writer.writerow(e)

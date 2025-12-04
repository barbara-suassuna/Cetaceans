import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
 
# Input/Output
input_fasta = "supermatrix_data1.fasta"
partitions_csv = "supermatrix_data1_partitions.csv"
output_nexus = "output_partitioned.nex"
 
# Load sequences
records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    record.annotations["molecule_type"] = "DNA"  # or RNA/protein
    records.append(record)
 
# Build NEXUS string manually to add partitions
nexus_lines = ["#NEXUS\n", "BEGIN DATA;"]
nexus_lines.append(f"    DIMENSIONS NTAX={len(records)} NCHAR={len(records[0].seq)};")
nexus_lines.append("    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX")
 
# Add sequences
for r in records:
    nexus_lines.append(f"    {r.id}    {r.seq}")
nexus_lines.append("    ;\nEND;\n")
 
# Read partitions and add a CHARSET block
nexus_lines.append("BEGIN SETS;")
with open(partitions_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row["part"]
        start = int(row["start"])
        end = int(row["stop"])
        nexus_lines.append(f"    CHARSET {name} = {start}-{end};")
nexus_lines.append("END;")
 
# Write NEXUS file
with open(output_nexus, "w") as out:
    out.write("\n".join(nexus_lines))
 
print("Partitioned NEXUS file created:", output_nexus)
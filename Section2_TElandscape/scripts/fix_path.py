path = "/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts/validate_mycomobilome.py"
with open(path, 'r') as f:
    lines = f.readlines()

with open(path, 'w') as f:
    for line in lines:
        if 'gff_path =' in line and 'filteredRepeats.gff' in line:
            f.write('        gff_path = info["dir"] / f"{info[\'prefix\']}_summaryFiles" / f"{info[\'prefix\']}.filteredRepeats.gff"\n')
        else:
            f.write(line)

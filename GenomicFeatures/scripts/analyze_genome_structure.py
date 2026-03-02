#!/usr/bin/env python3
"""
Análisis preciso de estructura genómica v4.0
CORRECCIONES:
- Maneja múltiples isoformas correctamente
- Usa features "gene" cuando disponibles (GFF3)
- Agrupa transcritos por gen correctamente
- CUENTA EXONES SOLO DE LA ISOFORMA MÁS LARGA (CORREGIDO)
- Reporta número de isoformas por gen
- NUEVO: Separa genes por biotype (protein_coding vs otros)
"""

import os
import sys
from collections import defaultdict
from Bio import SeqIO
import gzip
import re

def read_fasta_lengths(fasta_file):
    """Lee tamaños de scaffolds del FASTA"""
    scaffold_sizes = {}
    
    try:
        if fasta_file.endswith('.gz'):
            handle = gzip.open(fasta_file, 'rt')
        else:
            handle = open(fasta_file, 'r')
        
        for record in SeqIO.parse(handle, 'fasta'):
            scaffold_sizes[record.id] = len(record.seq)
        
        handle.close()
        return scaffold_sizes
    except Exception as e:
        print(f"    ERROR leyendo FASTA: {e}")
        return {}

def merge_intervals(intervals):
    """Fusiona intervalos solapantes"""
    if not intervals:
        return 0, []
    
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = []
    
    current_start, current_end = sorted_intervals[0]
    
    for start, end in sorted_intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    
    merged.append((current_start, current_end))
    
    total_length = sum(end - start + 1 for start, end in merged)
    return total_length, merged

def parse_attributes(attr_string):
    """Parse atributos GFF2/GFF3"""
    attrs = {}
    
    # GFF3 style: Key=Value
    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key.strip()] = value.strip()
    
    # GFF2 style: Key "Value"
    for match in re.finditer(r'(\w+)\s+"([^"]+)"', attr_string):
        attrs[match.group(1)] = match.group(2)
    
    # GFF2 style: Key Value (sin comillas)
    for match in re.finditer(r'(\w+)\s+([^\s;]+)', attr_string):
        if match.group(1) not in attrs:
            attrs[match.group(1)] = match.group(2)
    
    return attrs

def analyze_species(gff_file, fasta_file, species_name, output_dir):
    """Análisis principal de una especie"""
    
    print(f"\n{'='*60}")
    print(f"Procesando: {species_name}")
    print(f"{'='*60}")
    
    # Leer tamaños de scaffolds
    scaffold_sizes = read_fasta_lengths(fasta_file)
    
    if not scaffold_sizes:
        print(f"  ⚠️  No se pudo leer el FASTA")
        return
    
    total_genome_size = sum(scaffold_sizes.values())
    print(f"  📊 Genoma: {total_genome_size:,} bp en {len(scaffold_sizes)} scaffolds")
    
    # Estructuras de datos
    genes_from_feature = {}  # gene_id -> {'biotype', 'coords'}
    transcripts = {}  # transcript_id -> {'gene_id', 'biotype', 'scaffold', 'strand', 'exons': []}
    parent_mapping = {}  # mRNA/transcript -> gene
    
    # PRIMERA PASADA: Leer genes y relaciones
    try:
        if gff_file.endswith('.gz'):
            f = gzip.open(gff_file, 'rt')
        else:
            f = open(gff_file, 'r')
        
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqid, source, feat_type, start, end = fields[0], fields[1], fields[2], int(fields[3]), int(fields[4])
            strand, attributes = fields[6], fields[8]
            
            attrs = parse_attributes(attributes)
            
            # Capturar features "gene" directamente
            if feat_type == 'gene':
                gene_id = attrs.get('ID', attrs.get('gene_id', attrs.get('Name')))
                biotype = attrs.get('gene_biotype', attrs.get('biotype', 'unknown'))
                
                if gene_id:
                    genes_from_feature[gene_id] = {
                        'scaffold': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'biotype': biotype
                    }
            
            # Capturar mRNA/transcript y su Parent
            elif feat_type in ['mRNA', 'transcript']:
                transcript_id = attrs.get('ID', attrs.get('transcript_id', attrs.get('transcriptId')))
                parent_gene = attrs.get('Parent', attrs.get('gene_id'))
                biotype = attrs.get('transcript_biotype', attrs.get('biotype', 'unknown'))
                
                if transcript_id:
                    transcripts[transcript_id] = {
                        'gene_id': parent_gene,
                        'biotype': biotype,
                        'scaffold': seqid,
                        'strand': strand,
                        'exons': []
                    }
                    if parent_gene:
                        parent_mapping[transcript_id] = parent_gene
        
        f.close()
        
    except Exception as e:
        print(f"  ❌ ERROR en primera pasada: {e}")
        return
    
    # SEGUNDA PASADA: Leer exones
    try:
        if gff_file.endswith('.gz'):
            f = gzip.open(gff_file, 'rt')
        else:
            f = open(gff_file, 'r')
        
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqid, feat_type, start, end = fields[0], fields[2], int(fields[3]), int(fields[4])
            strand, attributes = fields[6], fields[8]
            
            # Solo usar 'exon', NO 'CDS' (para evitar duplicados en JGI)
            if feat_type != 'exon':
                continue
            
            attrs = parse_attributes(attributes)
            
            # Para JGI usar 'name' (gene/protein ID) como agrupador principal
            # Para otros usar Parent/transcript_id
            gene_id = attrs.get('name')  # JGI usa "name" para el gen/proteína
            transcript_id = attrs.get('transcriptId', attrs.get('transcript_id', attrs.get('Parent')))
            
            # Si no hay 'name' (no es JGI), usar Parent hierarchy
            if not gene_id:
                parent = attrs.get('Parent', attrs.get('gene_id', transcript_id))
                gene_id = parent
                transcript_id = parent
            
            if not gene_id:
                continue
            
            # Si el parent es un transcript conocido
            if transcript_id and transcript_id in transcripts:
                transcripts[transcript_id]['exons'].append((start, end))
            
            # Crear o actualizar usando gene_id como clave
            else:
                # Usar un composite key: (scaffold, gene_id) para evitar colisiones entre scaffolds
                composite_key = f"{seqid}:{gene_id}"
                
                if composite_key not in transcripts:
                    transcripts[composite_key] = {
                        'gene_id': gene_id,
                        'biotype': 'unknown',  # JGI no tiene biotypes
                        'scaffold': seqid,
                        'strand': strand,
                        'exons': []
                    }
                transcripts[composite_key]['exons'].append((start, end))
        
        f.close()
        
    except Exception as e:
        print(f"  ❌ ERROR en segunda pasada: {e}")
        return
    
    # AGRUPAR TRANSCRITOS POR GEN
    genes = defaultdict(lambda: {
        'scaffold': None,
        'strand': None,
        'biotype': 'unknown',
        'transcripts': [],
        'gene_coords': None
    })
    
    for transcript_id, t_data in transcripts.items():
        gene_id = t_data['gene_id'] if t_data['gene_id'] else transcript_id
        
        genes[gene_id]['scaffold'] = t_data['scaffold']
        genes[gene_id]['strand'] = t_data['strand']
        genes[gene_id]['biotype'] = t_data.get('biotype', 'unknown')
        genes[gene_id]['transcripts'].append(t_data['exons'])
        
        # Si tenemos coordenadas del gen desde feature "gene"
        if gene_id in genes_from_feature:
            genes[gene_id]['gene_coords'] = (
                genes_from_feature[gene_id]['start'],
                genes_from_feature[gene_id]['end']
            )
            # Usar biotype del gene si está disponible
            genes[gene_id]['biotype'] = genes_from_feature[gene_id].get('biotype', 'unknown')
    
    print(f"  🧬 Genes encontrados: {len(genes)}")
    print(f"  📝 Transcritos encontrados: {len(transcripts)}")
    
    # Contar biotypes
    biotype_counts = defaultdict(int)
    for gene_data in genes.values():
        biotype_counts[gene_data['biotype']] += 1
    
    print(f"  🏷️  Biotypes detectados:")
    for biotype, count in sorted(biotype_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"     - {biotype}: {count}")
    
    # ANÁLISIS POR GEN - SEPARADO POR BIOTYPE
    def analyze_gene_set(gene_subset, label):
        """Analiza un subconjunto de genes"""
        
        if not gene_subset:
            return None
        
        intron_counts = []
        intron_lengths = []
        exon_counts = []
        isoform_counts = []
        
        all_exon_intervals = defaultdict(list)
        all_intron_intervals = defaultdict(list)
        all_gene_intervals = defaultdict(list)
        
        genes_with_introns = 0
        total_introns = 0
        total_exons = 0
        total_isoforms = 0
        
        for gene_id, gene_data in gene_subset.items():
            scaffold = gene_data['scaffold']
            
            if not gene_data['transcripts']:
                continue
            
            # Contar número de isoformas
            num_isoforms = len(gene_data['transcripts'])
            total_isoforms += num_isoforms
            isoform_counts.append((gene_id, num_isoforms))
            
            # Elegir SOLO la isoforma más larga (mayor longitud en bp)
            longest_transcript = max(gene_data['transcripts'], key=lambda t: sum(e[1]-e[0]+1 for e in t))
            
            if not longest_transcript:
                continue
            
            # Ordenar exones
            exons_sorted = sorted(longest_transcript, key=lambda x: x[0])
            
            # Contar exones SOLO de la isoforma más larga
            num_exons = len(exons_sorted)
            total_exons += num_exons
            exon_counts.append((gene_id, num_exons))
            
            # Límites del gen
            if gene_data['gene_coords']:
                gene_start, gene_end = gene_data['gene_coords']
            else:
                gene_start = exons_sorted[0][0]
                gene_end = exons_sorted[-1][1]
            
            all_gene_intervals[scaffold].append((gene_start, gene_end))
            
            # Añadir exones
            for exon in exons_sorted:
                all_exon_intervals[scaffold].append(exon)
            
            # Calcular intrones
            num_introns = 0
            if len(exons_sorted) > 1:
                for i in range(len(exons_sorted) - 1):
                    intron_start = exons_sorted[i][1] + 1
                    intron_end = exons_sorted[i + 1][0] - 1
                    
                    if intron_end >= intron_start:
                        intron_length = intron_end - intron_start + 1
                        num_introns += 1
                        intron_lengths.append(intron_length)
                        all_intron_intervals[scaffold].append((intron_start, intron_end))
            
            if num_introns > 0:
                genes_with_introns += 1
            
            total_introns += num_introns
            intron_counts.append((gene_id, num_introns))
        
        # Calcular coberturas
        total_genic_bp = 0
        for scaffold in all_gene_intervals:
            bp, _ = merge_intervals(all_gene_intervals[scaffold])
            total_genic_bp += bp
        
        total_exon_bp = 0
        for scaffold in all_exon_intervals:
            bp, _ = merge_intervals(all_exon_intervals[scaffold])
            total_exon_bp += bp
        
        total_intron_bp = 0
        for scaffold in all_intron_intervals:
            bp, _ = merge_intervals(all_intron_intervals[scaffold])
            total_intron_bp += bp
        
        total_intergenic_bp = total_genome_size - total_genic_bp
        
        # Porcentajes
        pct_genic = (total_genic_bp / total_genome_size) * 100
        pct_exonic = (total_exon_bp / total_genome_size) * 100
        pct_intronic = (total_intron_bp / total_genome_size) * 100
        pct_intergenic = (total_intergenic_bp / total_genome_size) * 100
        
        # Estadísticas
        num_genes = len(gene_subset)
        avg_introns_all = total_introns / num_genes if num_genes > 0 else 0
        avg_introns_multiexon = total_introns / genes_with_introns if genes_with_introns > 0 else 0
        avg_exons_all = total_exons / num_genes if num_genes > 0 else 0
        avg_isoforms_per_gene = total_isoforms / num_genes if num_genes > 0 else 0
        
        return {
            'label': label,
            'num_genes': num_genes,
            'total_isoforms': total_isoforms,
            'total_exons': total_exons,
            'total_introns': total_introns,
            'genes_with_introns': genes_with_introns,
            'avg_isoforms': avg_isoforms_per_gene,
            'avg_exons': avg_exons_all,
            'avg_introns_all': avg_introns_all,
            'avg_introns_multiexon': avg_introns_multiexon,
            'pct_genic': pct_genic,
            'pct_exonic': pct_exonic,
            'pct_intronic': pct_intronic,
            'pct_intergenic': pct_intergenic,
            'intron_counts': intron_counts,
            'exon_counts': exon_counts,
            'isoform_counts': isoform_counts,
            'intron_lengths': intron_lengths
        }
    
    # Analizar todos los genes
    all_results = analyze_gene_set(genes, 'ALL')
    
    # Analizar protein_coding genes
    protein_coding_genes = {gid: gdata for gid, gdata in genes.items() 
                           if gdata['biotype'] == 'protein_coding'}
    pc_results = analyze_gene_set(protein_coding_genes, 'PROTEIN_CODING')
    
    # Analizar non-coding genes
    non_coding_genes = {gid: gdata for gid, gdata in genes.items() 
                       if gdata['biotype'] != 'protein_coding' and gdata['biotype'] != 'unknown'}
    nc_results = analyze_gene_set(non_coding_genes, 'NON_CODING')
    
    print(f"  🔬 Calculando coberturas genómicas...")
    
    # Mostrar resultados
    for results in [all_results, pc_results, nc_results]:
        if results:
            print(f"\n  📈 RESULTADOS [{results['label']}]:")
            print(f"     Genes totales: {results['num_genes']:,}")
            print(f"     Isoformas totales: {results['total_isoforms']:,}")
            print(f"     Isoformas/gen (promedio): {results['avg_isoforms']:.2f}")
            print(f"     Exones totales (isoforma larga): {results['total_exons']:,}")
            print(f"     Exones/gen (promedio): {results['avg_exons']:.2f}")
            print(f"     Intrones totales: {results['total_introns']:,}")
            print(f"     Intrones/gen (todos): {results['avg_introns_all']:.2f}")
            print(f"     Intrones/gen (multiexón): {results['avg_introns_multiexon']:.2f}")
    
    # Mostrar cobertura solo para ALL
    if all_results:
        print(f"\n  🎯 COBERTURA GENÓMICA [ALL]:")
        print(f"     Génica:       {all_results['pct_genic']:6.2f}%")
        print(f"     - Exones:     {all_results['pct_exonic']:6.2f}%")
        print(f"     - Intrones:   {all_results['pct_intronic']:6.2f}%")
        print(f"     Intergénica:  {all_results['pct_intergenic']:6.2f}%")
    
    # Guardar resultados
    os.makedirs(output_dir, exist_ok=True)
    
    # Guardar para cada categoría
    for results in [all_results, pc_results, nc_results]:
        if results:
            suffix = f"_{results['label']}"
            
            # 1. Resumen
            with open(os.path.join(output_dir, f'{species_name}{suffix}_SUMMARY.txt'), 'w') as f:
                f.write("Metric\tValue\n")
                f.write(f"Genome_Size_bp\t{total_genome_size}\n")
                f.write(f"Total_Genes\t{results['num_genes']}\n")
                f.write(f"Total_Isoforms\t{results['total_isoforms']}\n")
                f.write(f"Avg_Isoforms_Per_Gene\t{results['avg_isoforms']:.2f}\n")
                f.write(f"Total_Exons\t{results['total_exons']}\n")
                f.write(f"Total_Introns\t{results['total_introns']}\n")
                f.write(f"Avg_Exons_Per_Gene\t{results['avg_exons']:.2f}\n")
                f.write(f"Avg_Introns_All\t{results['avg_introns_all']:.2f}\n")
                f.write(f"Avg_Introns_Multiexon\t{results['avg_introns_multiexon']:.2f}\n")
                f.write(f"Pct_Genic\t{results['pct_genic']:.4f}\n")
                f.write(f"Pct_Exonic\t{results['pct_exonic']:.4f}\n")
                f.write(f"Pct_Intronic\t{results['pct_intronic']:.4f}\n")
                f.write(f"Pct_Intergenic\t{results['pct_intergenic']:.4f}\n")
            
            # 2. Intrones por gen
            with open(os.path.join(output_dir, f'{species_name}{suffix}_intron_counts.txt'), 'w') as f:
                f.write("Gene\tIntron_count\n")
                for gene_id, count in results['intron_counts']:
                    f.write(f"{gene_id}\t{count}\n")
            
            # 3. Exones por gen
            with open(os.path.join(output_dir, f'{species_name}{suffix}_exon_counts.txt'), 'w') as f:
                f.write("Gene\tExon_count\n")
                for gene_id, count in results['exon_counts']:
                    f.write(f"{gene_id}\t{count}\n")
            
            # 4. Longitudes de intrones
            with open(os.path.join(output_dir, f'{species_name}{suffix}_intron_sizes.txt'), 'w') as f:
                f.write("Intron_length\n")
                for length in results['intron_lengths']:
                    f.write(f"{length}\n")
            
            # 5. Isoformas por gen
            with open(os.path.join(output_dir, f'{species_name}{suffix}_isoform_counts.txt'), 'w') as f:
                f.write("Gene\tIsoform_count\n")
                for gene_id, count in results['isoform_counts']:
                    f.write(f"{gene_id}\t{count}\n")
    
    print(f"\n  💾 Archivos guardados en: {output_dir}")


def main():
    if len(sys.argv) < 4:
        print("Uso: python script.py <dir_gff> <dir_genomas> <dir_resultados>")
        sys.exit(1)
    
    gff_dir = sys.argv[1]
    genome_dir = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Encontrar archivos GFF
    import glob
    gff_files = []
    for ext in ['*.gff', '*.gff3', '*.gff.gz', '*.gff3.gz']:
        gff_files.extend(glob.glob(os.path.join(gff_dir, ext)))
    
    print(f"\n🔍 Encontrados {len(gff_files)} archivos GFF")
    
    # Procesar cada especie
    for gff_file in sorted(gff_files):
        base_name = os.path.basename(gff_file)
        species = base_name.replace('.gff3', '').replace('.gff', '').replace('.gz', '')
        
        # Buscar FASTA correspondiente
        fasta_file = None
        for ext in ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz']:
            candidate = os.path.join(genome_dir, species + ext)
            if os.path.exists(candidate):
                fasta_file = candidate
                break
        
        if not fasta_file:
            print(f"\n⚠️  No se encontró FASTA para {species}")
            continue
        
        analyze_species(gff_file, fasta_file, species, output_dir)
    
    # Consolidar summaries - PARA CADA CATEGORÍA
    print(f"\n{'='*60}")
    print("📊 Consolidando resultados...")
    print(f"{'='*60}")
    
    for category in ['_ALL', '_PROTEIN_CODING', '_NON_CODING']:
        summary_files = glob.glob(os.path.join(output_dir, f'*{category}_SUMMARY.txt'))
        
        if not summary_files:
            continue
        
        consolidated_data = []
        for summary_file in sorted(summary_files):
            species = os.path.basename(summary_file).replace(f'{category}_SUMMARY.txt', '')
            
            data = {'Species': species}
            with open(summary_file, 'r') as f:
                f.readline()  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        metric, value = parts
                        data[metric] = value
            
            consolidated_data.append(data)
        
        # Escribir consolidado
        output_file = f'ALL_SPECIES{category}_SUMMARY.txt'
        with open(os.path.join(output_dir, output_file), 'w') as f:
            headers = ['Species', 'Genome_Size_bp', 'Total_Genes', 'Total_Isoforms', 'Avg_Isoforms_Per_Gene',
                      'Total_Exons', 'Total_Introns', 'Avg_Exons_Per_Gene', 'Avg_Introns_All', 
                      'Avg_Introns_Multiexon', 'Pct_Genic', 'Pct_Exonic', 'Pct_Intronic', 'Pct_Intergenic']
            f.write('\t'.join(headers) + '\n')
            
            for data in consolidated_data:
                row = [data.get(h, 'NA') for h in headers]
                f.write('\t'.join(str(x) for x in row) + '\n')
        
        print(f"   ✅ {output_file}")
    
    print(f"\n✅ Análisis completo!")

if __name__ == '__main__':
    main()
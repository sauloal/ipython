#!/usr/bin/python

import os
import sys
from sqlalchemy.util import KeyedTuple

#grep -P '#|\tgap' SL2.50ch_from_sc.agp.gff3 > SL2.50ch_from_sc.agp.gff3.gap.gff3
#grep -P '#|\tcontig' SL2.50ch_from_sc.agp.gff3 > SL2.50ch_from_sc.agp.gff3.contig.gff3
#grep -P '#|\tscaffold' SL2.50ch_from_sc.agp.gff3 > SL2.50ch_from_sc.agp.gff3.scaffold.gff3


flags_gap = ('N', 'U')
col_names_shared =                    [ 'object_id'   , 'object_beg'   , 'object_end'   , 'part_number'     , 'component_type' ]
col_names_gap_Y  = col_names_shared + [ 'gap_length'  , 'gap_type'     , 'linkage'      , 'linkage_evidence'                   ]
col_names_gap_N  = col_names_shared + [ 'gap_length'  , 'gap_type'     , 'linkage'                                             ]
col_names_contig = col_names_shared + [ 'component_id', 'component_beg', 'component_end', 'orientation'                        ]

component_types = {
    'A': 'Active Finishing',
    'D': 'Draft HTG (often phase1 and phase2 are called Draft, whether or not they have the draft keyword)',
    'F': 'Finished HTG (phase3)',
    'G': 'Whole Genome Finishing',
    'O': 'Other sequence (typically means no HTG keyword)',
    'P': 'Pre Draft',
    'W': 'WGS contig',
    'N': 'gap with specified size',
    'U': 'gap of unknown size, defaulting to 100 bases'
}

type_mapper = {
    'A': 'contig',
    'D': 'contig',
    'F': 'scaffold',
    'G': 'scaffold',
    'O': 'contig',
    'P': 'contig',
    'W': 'contig',
    'N': 'gap',
    'U': 'gap'
}

gff_cols    = [ 'object_id', 'source_name', 'type', 'object_beg', 'object_end', 'score', 'orientation', 'phase' ]
source_name = 'agp'
id_prefix   = 'agp_'

def main(args):
    in_agp = args[0]
    ou_gff = in_agp + ".gff3"
    
    print "SAVING TO", ou_gff
    with open(in_agp, 'r') as fhd_in:
        with open(ou_gff, 'w') as fhd_ou:
            fhd_ou.write("##gff-version 3\n")
            fhd_ou.write("#infile: %s\n" % in_agp)
            
            for line_in in fhd_in:
                line_in = line_in.strip()
                
                if len(line_in) == 0:
                    continue
                
                if line_in[0] == "#":
                    continue
                
                cols_in = line_in.split("\t")
                
                if cols_in[4] in flags_gap:
                    if   cols_in[7] == 'no':
                        assert len(cols_in) == len(col_names_gap_N), "%d vs %d: %s" % ( len(cols_in), len(col_names_gap_N), " ".join(cols_in) )
                        cols_in = KeyedTuple(cols_in, labels=col_names_gap_N   )._asdict()
                    
                    elif cols_in[7] == 'yes':
                        assert len(cols_in) == len(col_names_gap_Y), "%d vs %d: %s" % ( len(cols_in), len(col_names_gap_Y), " ".join(cols_in) )
                        cols_in = KeyedTuple(cols_in, labels=col_names_gap_Y   )._asdict()

                    else:
                        print "unknown gap linkage evidence:", cols_in[7]
                        sys.exit(1)

                    cols_in['orientation'        ] = '.'
                
                else:
                    assert len(cols_in) == len(col_names_contig), "%d vs %d: %s" % ( len(cols_in), len(col_names_contig), " ".join(cols_in) )
                    cols_in = KeyedTuple(cols_in, labels=col_names_contig)._asdict()
                
                cols_in['component_type_desc'] = component_types[cols_in['component_type']]
                cols_in['type'               ] = type_mapper[    cols_in['component_type']]
                cols_in['source_name'        ] = source_name
                cols_in['score'              ] = '.'
                cols_in['phase'              ] = '.'
                
                
                #print cols_in
                row_id = id_prefix + cols_in['object_id'] + '_' + cols_in['part_number']
                
                fhd_ou.write("\t".join( [cols_in[x] for x in gff_cols] ) )
                fhd_ou.write("\tID=%s;Name=%s;" % ( row_id, row_id ) )
                fhd_ou.write(";".join( ["=".join([x, cols_in[x]]) for x in sorted(cols_in) if x not in gff_cols] ) )
                fhd_ou.write("\n")


if __name__ == '__main__':
    main(sys.argv[1:])
        
#1               2        3               4       5       6               7       8        9
#SL2.50ch00      1        2191949         1       W       SL2.40sc05082   1       2191949  0
#SL2.50ch00      2191950  2192049         2       U       100             contig  no
#SL2.50ch00      21801621 21801720        6260    U       100             contig  no
#SL2.50ch00      21801721 21803721        6261    W       SL2.40sc03931   1       2001     0
#SL2.50ch00      21803722 21803821        6262    U       100             contig  no
#SL2.50ch00      21803822 21805821        6263    W       SL2.40sc04627   1       2000     0
#SL2.50ch01      1        32987597        1       F       SL2.40sc04133   1       32987597        -
#SL2.50ch01      32987598 35267597        2       N       2280000         contig  no
#SL2.50ch01      35267598 36990194        3       F       SL2.40sc04191   1       1722597  +



#http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml


#1 object This is the identifier for the object being assembled. This can be a
#chromosome, scaffold or contig. If an accession.version identifier is not used
#to describe the object the naming convention is to precede chromosome numbers
#with 'chr' (e.g. chr1) and linkage group numbers with 'LG' (e.g. LG3). Contigs
#or scaffolds may have any identifier that is unique within the assembly.

#2 object_beg The starting coordinates of the component/gap on the object in
#column 1. These are the location in the object's coordinate system, not the
#component's.

#3 object_end The ending coordinates of the component/gap on the object in
#column 1. These are the location in the object's coordinate system, not the
#component's.

#4 part_number The line count for the components/gaps that make up the object
#described in column 1.

#5 component_type The sequencing status of the component. These typically
#correspond to keywords in the International Sequence Database
#(GenBank/EMBL/DDBJ) submission. Current acceptable values are:

#The sequencing status of the component. These typically correspond to keywords
#in the International Sequence Database (GenBank/EMBL/DDBJ) submission. Current
#acceptable values are:
#A        Active Finishing
#D        Draft HTG (often phase1 and phase2 are called Draft, whether or not they have the draft keyword).
#F        Finished HTG (phase3)
#G        Whole Genome Finishing
#O        Other sequence (typically means no HTG keyword)
#P        Pre Draft
#W        WGS contig
#N        gap with specified size
#U        gap of unknown size, defaulting to 100 bases.

#6a component_id If column 5 not equal to N or U: This is a unique identifier
#for the sequence component contributing to the object described in column 1.
#Ideally this will be a valid accession.version identifier as assigned by
#GenBank/EMBL/DDBJ. If the sequence has not been submitted to a public
#repository yet a local identifier should be used.

#6b gap_length If column 5 equal to N or U: This column represents the length of
#the gap. N type gaps can be of any length. A length of 100 must be used for all
#U type gaps.

#7a component_beg If column 5 not equal to N or U: This column specifies the
#beginning of the part of the component sequence that contributes to the object
#in column 1 (in component coordinates).

#7b	gap_type	
#    If column 5 equal to N or U: This column specifies the gap type.
#    Accepted values:
#        scaffold: a gap between two sequence contigs in a scaffold (superscaffold or ultra-scaffold).
#        contig: an unspanned gap between two sequence contigs.
#        centromere: a gap inserted for the centromere.
#        short_arm: a gap inserted at the start of an acrocentric chromosome.
#        heterochromatin: a gap inserted for an especially large region of heterochromatic sequence (may also include the centromere).
#        telomere: a gap inserted for the telomere.
#        repeat: an unresolvable repeat.

#8a component_end If column 5 not equal to N or U: This column specifies the end
#of the part of the component that contributes to the object in column 1 (in
#component coordinates).

#8b	linkage	        If column 5 equal to N or U: This column indicates if there is evidence of linkage between the adjacent lines.
#    Values: yes, no

#9a	orientation	
#    If column 5 not equal to N or U: This column specifies the orientation of the component relative to the object in column 1.
#    Values:
#        + plus
#        - minus
#        ? unknown
#        0 (zero) unknown (deprecated)
#        na irrelevant
#    By default, components with unknown orientation (?, 0 or na) are treated as if they had + orientation.

#9b Linkage evidence If column 5 equal to N or U: This specifies the type of
#evidence used to assert linkage (as indicated in column 8b).
#    Accepted values:
#        na            used when no linkage is being asserted (column 8b is 'no')
#        paired-ends   paired sequences from the two ends of a DNA fragment.
#        align_genus   alignment to a reference genome within the same genus.
#        align_xgenus  alignment to a reference genome within another genus.
#        align_trnscpt alignment to a transcript from the same species.
#        within_clone sequence on both sides of the gap is derived from the same
#                     clone, but the gap is not spanned by paired-ends. The adjacent sequence
#                     contigs have unknown order and orientation.
#        clone_contig linkage is provided by a clone contig in the tiling path
#                     (TPF). For example, a gap where there is a known clone, but there is
#                     not yet sequence for that clone.
#        map           linkage asserted using a non-sequence based map such as RH, linkage, fingerprint or optical.
#        strobe        strobe sequencing (PacBio).
#        unspecified   used only when converting old AGPs that lack a field for linkage evidence into the new format.
#    If there are multiple lines of evidence to support linkage, all can be
#    listed using a ';' delimiter (e.g. paired-ends;align_xgenus).



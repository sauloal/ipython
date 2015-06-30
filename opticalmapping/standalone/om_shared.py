import re
import operator
import argparse
import textwrap
from collections     import defaultdict
from sqlalchemy.util import KeyedTuple

"""
#h XmapEntryID  QryContigID     RefContigID     QryStartPos     QryEndPos       RefStartPos     RefEndPos       Orientation     Confidence      HitEnum QryLen     RefLen    LabelChannel    Alignment
#f int          int             int             float           float           float           float           string          float           string  float      float     int             string
   1            141             1               528400.6        571697.5        10672           54237.5         +               6.65            4M2D2M  1439123.5  21805821  1               "(1,34)(2,34)(3,35)(4,36)(5,37)(6,38)(8,38)(9,39)"
"""

cols_to_index  = ["QryContigID", "RefContigID", "Orientation", "XmapEntryID"]
group_by       = [
    ["RefContigID", "QryContigID"],
    ["RefContigID", "RefStartPos"],
    ["RefContigID", "RefEndPos"  ],
    ["QryContigID", "RefContigID"],
    ["QryContigID", "Orientation"],
    ["QryContigID", "XmapEntryID"],
    ["XmapEntryID", "Confidence" ],
]

re_matches    = re.compile("(\d+)M")
re_insertions = re.compile("(\d+)I")
re_deletions  = re.compile("(\d+)D")
re_alignment  = re.compile("\((\d+),(\d+)\)")



def col_parse_orientation(val):
    assert(val in ("+", "-"))
    return val

def col_parse_hit_enum(val):
    #4M2D2M
    assert(set([x for x in val]) <= set(['M', 'D', 'I', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']))
    return val

def col_parse_alignment(val):
    #"(1,34)(2,34)(3,35)(4,36)(5,37)(6,38)(8,38)(9,39)"
    val = val.strip('"')
    assert(set([x for x in val]) <= set(['(', ')', ',', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']))
    return val

def col_parse_bool(val):
    if   val.lower() in ("t", "true" , "1" ):
        return True
    
    elif val.lower() in ("f", "false", "0"):
        return False
    
    else:
        print "could not parse bool %s" % val
        sys.exit(1)





def process_cigar(cigar):
    """
    2M3D1M1D1M1D4M1I2M1D2M1D1M2I2D9M3I3M1D6M1D2M2D1M1D6M1D1M1D1M2D2M2D1M1I1D1M1D5M2D4M2D1M2D2M1D2M1D3M1D1M1D2M3I3D1M1D1M3D2M3D1M2I1D1M2D1M1D1M1I2D3M2I1M1D2M1D1M1D1M2I3D3M3D1M2D1M1D1M1D5M2D12M
    """
    assert(set([x for x in cigar]) <= set(['M', 'D', 'I', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']))

    cigar_matches    = 0
    cigar_insertions = 0
    cigar_deletions  = 0

    i_matches = re_matches   .finditer(cigar)
    i_inserts = re_insertions.finditer(cigar)
    i_deletes = re_deletions .finditer(cigar)

    for i in i_matches:
        n                 = i.group(1)
        cigar_matches    += int(n)

    for i in i_inserts:
        n                 = i.group(1)
        cigar_insertions += int(n)

    for i in i_deletes:
        n                 = i.group(1)
        cigar_deletions  += int(n)
    
    return cigar_matches, cigar_insertions, cigar_deletions

def process_alignment(alignment):
    """
     Alignment (4862,48)(4863,48)(4864,47)(4865,46)(4866,45)(4867,44)(4870,43)(4873,42)(4874,41)(4875,40)(4877,40)(4878,39)(4879,38)(4880,37)(4883,36)(4884,36)(4885,35)(4886,34)(4887,33)(4888,33)(4889,32)(4890,30)(4891,30)(4892,29)(4893,28)(4894,28)(4899,27)(4900,26)(4901,25)(4902,24)(4903,23)(4904,22)(4906,21)(4907,21)(4908,20)(4910,19)(4911,18)(4912,17)(4913,16)(4915,15)(4917,14)(4918,13)(4919,12)(4920,11)(4922,10)(4923,9)(4925,8)(4927,7)(4930,6)(4931,5)(4932,3)(4933,2)(4934,1)
    """
    assert(set([x for x in alignment]) <= set(['(', ')', ',', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']))
    
    count_refs             = defaultdict(int)
    count_queries          = defaultdict(int)
    count_refs_colapses    = 0
    count_queries_colapses = 0
    
    i_alignment = re_alignment.finditer(alignment)
    for i in i_alignment:
        c_r           = int(i.group(1))
        c_q           = int(i.group(2))

        count_refs   [c_r] += 1
        count_queries[c_q] += 1

    count_refs_colapses    = sum([count_refs[   x] for x in count_refs    if count_refs[   x] > 1])
    count_queries_colapses = sum([count_queries[x] for x in count_queries if count_queries[x] > 1])

    return len(count_refs), len(count_queries), count_refs_colapses, count_queries_colapses




def gen_valid_fields(valid_fields):
    valid_fields['names'  ] = [None] * len(valid_fields['data'])
    valid_fields['parsers'] = {}
    valid_fields['types'  ] = {}
    valid_fields['poses'  ] = {}
    valid_fields['helps'  ] = {}
    valid_fields['helps_t'] = {}
    
    for field_pos, field_data, in enumerate(valid_fields['data']):
        field_name, field_type, field_parser, field_help = field_data
        
        valid_fields['names'  ][ field_pos  ] = field_name
        valid_fields['parsers'][ field_name ] = field_parser
        valid_fields['types'  ][ field_name ] = field_type
        valid_fields['poses'  ][ field_name ] = field_pos
        valid_fields['helps'  ][ field_name ] = ("\n"  + (" "*53)).join(textwrap.wrap(field_help, 80))
        valid_fields['helps_t'][ field_name ] = ("\n#" + (" "*42)).join(textwrap.wrap(field_help, 80))

    return valid_fields

def parser_in(orig_parser, value, sep=","):
    return set([ orig_parser(v) for v in value.split(sep) ])

def function_in(val, qry):
    return val in qry

valid_operators = {
    #op name      op func            override parser   help
    'eq'      : [ operator.eq      , None, "value <Equal> to filter" ],
    'ge'      : [ operator.ge      , None, "value <Greater than or Equal> to filter" ],
    'gt'      : [ operator.gt      , None, "value <Greater than> filter" ],
    #'is_not'  : [ operator.is_not  , None, "value <Is not> filter [class]" ],
    #'is'      : [ operator.is_     , None, "value <Is> filter [class]" ],
    'le'      : [ operator.le      , None, "value <Less than or Equal> to filter" ],
    'lt'      : [ operator.lt      , None, "value <Less than> filter" ],
    'ne'      : [ operator.ne      , None, "value <Not equal> to filter" ],
    #'truth'   : [ operator.truth   , None, "value is <Truth>" ],
    'contains': [ operator.contains, parser_in, "value <Contains> filter" ],
    'in'      : [ function_in      , parser_in, "value <In> filter [comman separated]" ]
}
#for operator_name in sorted(dir(operator)):
#    if operator_name[:2] == "__":
#        continue
#    if operator_name[-1] == "_":
#        continue
#    print "adding operator name %s" % operator_name
#    valid_operators[operator_name] = getattr(operator, operator_name)


def gen_filter(filter_datas, valid_fields):
    filters = []
    
    if filter_datas is not None:
        for filter_data in filter_datas:
            filter_cols = filter_data.split(":")

            if len(filter_cols) != 3:
                print "filter has to have 3 parts <field>:<function>:<value>, %d found in %s" % ( len(filter_cols), filter_data )
                sys.exit(0)

            field_name, operator_name, value_str = filter_cols
            assert field_name    in valid_fields['names'  ], "invalid value for field name"
            assert operator_name in valid_operators        , "operator %s does not exists. acceptable values are: lt, le, eq, ne, ge, gt" % operator_name

            if valid_operators[operator_name][1] is None:
                value          =                                    valid_fields['parsers'][field_name]( value_str )
            else:
                value          = valid_operators[operator_name][1]( valid_fields['parsers'][field_name], value_str )

            operator_val   = valid_operators[operator_name][0]
            filter_res     = [field_name, operator_name, operator_val, value_str, value]

            filters.append(filter_res)

    return filters



def parse_file(infile, valid_fields):
    data            = []
    names           = []
    seman           = {}
    types           = []
    indexer         = {}
    headers         = []
    filters         = []
    groups          = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    ref_maps_from   = ""
    query_maps_from = ""

    
    for cti in cols_to_index:
        indexer[cti] = defaultdict(set)
    
    with open(infile, 'r') as fhd:
        for line in fhd:
            line = line.strip()
            
            if len(line) == 0:
                continue
            
            if line[0] == "#":
                headers.append(line)
                
                
                if len(line) == 1:
                    continue
                
                if line[2:10] == "FILTER :":
                    print "PARSING FILTER"
                    #"Confidence                             :  ge : 10.0"
                    cols = [ x.strip() for x in line[2:].split(":") ]
                    filter_line = ":".join(cols[1:])
                    filters.append( filter_line )
                
                elif line[1] == "h":
                    line  = line[3:]
                    names = [x.strip() for x in line.split("\t")]
                    #print "NAMES", names
                    
                    for p, n in enumerate(names):
                        seman[n] = p
                
                elif line[1] == "f":
                    line  = line[3:]
                    types = [x.strip() for x in line.split("\t")]
                    
                    #for tp in xrange(len(types)):
                    #    t = types[tp]
                    #    if   t == "int":
                    #        types[tp] = int
                    #    elif t == "float":
                    #        types[tp] = float
                    #    elif t == "string":
                    #        types[tp] = col_parsers[ names[tp] ]
                    
                    assert(len(types) == len(names))
                    
                    #print "TYPES", types
                
                elif "Reference Maps From:" in line:
                    # Reference Maps From:  S_lycopersicum_chromosomes.2.50.BspQI-BbvCI_to_EXP_REFINEFINAL1_r.cmap
                    # Query Maps From:      S_lycopersicum_chromosomes.2.50.BspQI-BbvCI_to_EXP_REFINEFINAL1_q.cmap
                    ref_maps_from   = line[23:].strip()
    
                elif "Query Maps From:" in line:
                    query_maps_from = line[23:].strip()

                continue
            
            cols = [x.strip()                                  for x in line.split("\t") ]
            vals = [valid_fields['parsers'][names[p]](cols[p]) for p in xrange(len(cols))]

            for ind in indexer:
                indexer[ ind ][ vals[seman[ind]] ].add(len(data))
            
            for grp_from, grp_to in group_by:
                val_from = vals[seman[grp_from]]
                val_to   = vals[seman[grp_to  ]]
                groups[grp_from+'_'+grp_to][val_from][val_to].add(len(data))
            
            data.append(vals)

    return data, headers, names, seman, types, indexer, groups, ref_maps_from, query_maps_from, filters



def stats_from_data_vals(RefContigID, QryContigID, groups, indexer, data, data_vals, valid_data_poses):
    ref_lens             = [ ( x["RefStartPos"], x["RefEndPos"] ) for x in data_vals ]
    qry_lens             = [ ( x["QryStartPos"], x["QryEndPos"] ) for x in data_vals ]
    
    num_qry_matches = []
    for RefContigID_l in groups["QryContigID_RefContigID"][QryContigID]:
        for match_pos in groups["QryContigID_RefContigID"][QryContigID][RefContigID_l]:
            if match_pos in valid_data_poses:
                num_qry_matches.append(RefContigID_l)

    #num_qry_matches      = len( groups["QryContigID_RefContigID"][QryContigID] )
    num_qry_matches      = len( set(num_qry_matches)                       )
    num_orientations     = len( set([x["Orientation"] for x in data_vals]) )

    ref_no_gap_len       = sum( [ max(x)-min(x) for x in ref_lens ] )
    ref_min_coord        = min( [ min(x)        for x in ref_lens ] )
    ref_max_coord        = max( [ max(x)        for x in ref_lens ] )
    ref_gap_len          = ref_max_coord - ref_min_coord
    
    qry_no_gap_len       = sum( [ max(x)-min(x) for x in qry_lens ] )
    qry_min_coord        = min( [ min(x)        for x in qry_lens ] )
    qry_max_coord        = max( [ max(x)        for x in qry_lens ] )
    qry_gap_len          = qry_max_coord - qry_min_coord
    
    XmapEntryIDs         = groups["QryContigID_XmapEntryID"][QryContigID].keys()
    
    Confidences          = []
    for XmapEntryID in XmapEntryIDs:
        data_pos = list(indexer["XmapEntryID"][XmapEntryID])[0]
        if data_pos not in valid_data_poses:
            continue
        Confidences.append( [ data[data_pos]["Confidence"], data[data_pos]["RefContigID"] ] )
    
    max_confidence       = max([ x[0] for x in Confidences ])
    max_confidence_chrom = [ x[1] for x in Confidences if x[0] == max_confidence][0]

    stats = {}
    stats["_meta_is_max_confidence_for_qry_chrom"      ] = max_confidence_chrom == RefContigID
    
    stats["_meta_len_ref_match_gapped"                 ] = ref_gap_len
    stats["_meta_len_ref_match_no_gap"                 ] = ref_no_gap_len
    stats["_meta_len_qry_match_gapped"                 ] = qry_gap_len
    stats["_meta_len_qry_match_no_gap"                 ] = qry_no_gap_len

    stats["_meta_max_confidence_for_qry"               ] = max_confidence
    stats["_meta_max_confidence_for_qry_chrom"         ] = max_confidence_chrom

    stats["_meta_num_orientations"                     ] = num_orientations
    stats["_meta_num_qry_matches"                      ] = num_qry_matches
    stats["_meta_qry_matches"                          ] = ','.join( [ str(x) for x in sorted(list(set([ x[1] for x in Confidences ]))) ] )

    stats["_meta_proportion_sizes_gapped"              ] = (ref_gap_len    * 1.0)/ qry_gap_len
    stats["_meta_proportion_sizes_no_gap"              ] = (ref_no_gap_len * 1.0)/ qry_no_gap_len

    return stats

valid_fields_g  = {
    'data':
        [
            #colum name                                type        parser order help
            [ "XmapEntryID"                           , 'int'   , int                  , 'A unique line number for the data lines in the XMAP file. Note: For 2-color, the XmapEntryID will begin with the number 2.' ],
            
            [ "QryContigID"                           , 'int'   , int                  , 'Map ID of query map (Contig ID from .cmap file for query)' ],
            [ "RefContigID"                           , 'int'   , int                  , 'Map ID of the reference map from the .cmap reference file (the .cmap file may contain multiple reference maps). Note: RefContigIDs must be integers, but they need not be sequential.' ],
            
            [ "QryStartPos"                           , 'float' , float                , 'Coordinates of the first aligned label on the query map (Start position of hit on query map)' ],
            [ "QryEndPos"                             , 'float' , float                , 'Coordinates of the last aligned label on the query map (Stop position of hit on query map)' ],
            [ "RefStartPos"                           , 'float' , float                , 'Coordinates of the first aligned label on the reference or anchor map' ],
            [ "RefEndPos"                             , 'float' , float                , 'Coordinates of the last aligned label on the reference or anchor map' ],
        
            [ "Orientation"                           , 'string', col_parse_orientation, 'The relative orientation of the query map relative to the reference: forward (+) or reverse (-). The convention is that the reference is always positive orientation, so if the query aligns in reverse, it is shown as having negative (-) orientation. Note: For 2-color, the orientation will be the same.' ],
            [ "Confidence"                            , 'float' , float                , 'Statistical Confidence of result: Negative Log10 of p-value of alignment (without Bonferroni Correction for multiple experiments). Note: For 2-color, the confidence number is the combined confidence of the alignment for both colors.' ],
            [ "HitEnum"                               , 'string', col_parse_hit_enum   , 'Pseudo-CIGAR string representing matches (M), insertions (I), or deletions (D) of label sites with respect to the reference or anchor map. Count begins at the leftmost anchor label of that color. Note: When 2 or more anchor sites resolve into a single query site, only the rightmost anchor site is shown matched with the query site and the leftmost associated anchor sites are shown as deletions.' ],
            [ "QryLen"                                , 'float' , float                , 'Length of query map from _q.cmap.' ],
            [ "RefLen"                                , 'float' , float                , 'Length of reference map from _r.cmap.' ],
            [ "LabelChannel"                          , 'int'   , int                  , 'Color channel of alignment from cmap files. For 1-color data, LabelChannel is 1. For 2-color data: Using -usecolor N, the LabelChannel is N (N = 1 or 2), and there is only one XMAP entry per alignment for the color channel specified by N. Without -usecolor N, LabelChannel is 1 or 2. In this case, there are two XMAP entries (two lines), one for each color channel.' ],
            [ "Alignment"                             , 'string', col_parse_alignment  , 'Indices of the aligned site ID pairs. (When the query orientation is reversed ("-"), the query IDs are in descending order.) Count begins at the leftmost anchor label of that color. Note: When two sites in the reference align with the same site in the query, it is an indication that the two sites in the reference failed to resolve. Alignment provides a view of aligned pairs which would normally be ignored by HitEnum (CIGAR string).' ],
        
            [ "_meta_alignment_count_queries"         , 'int'   , int                  , 'Number of query labels in alignment' ],
            [ "_meta_alignment_count_queries_colapses", 'int'   , int                  , 'Number of query label collapses in alignment. A collapse happens when a label matches more than once a reference label' ],
            [ "_meta_alignment_count_refs"            , 'int'   , int                  , 'Number of reference labels in alignment' ],
            [ "_meta_alignment_count_refs_colapses"   , 'int'   , int                  , 'Number of reference label collapses in alignment. A collapse happens when a label matches more than once a query label' ],
            
            [ "_meta_cigar_deletions"                 , 'int'   , int                  , 'Number of deleted labels in CIGAR string' ],
            [ "_meta_cigar_insertions"                , 'int'   , int                  , 'Number of inserted labels in CIGAR string' ],
            [ "_meta_cigar_matches"                   , 'int'   , int                  , 'Number of match labels in CIGAR string' ],
            
            [ "_meta_is_max_confidence_for_qry_chrom" , 'bool'  , col_parse_bool       , 'Whether the current RefContigID is the highest confidence match for this QryContigID' ],
            
            [ "_meta_len_qry_match_gapped"            , 'float' , float                , 'Length of the query match including gaps' ],
            [ "_meta_len_qry_match_no_gap"            , 'float' , float                , 'Length of the query match excluding gaps' ],
            [ "_meta_len_ref_match_gapped"            , 'float' , float                , 'Length of the reference match including gaps' ],
            [ "_meta_len_ref_match_no_gap"            , 'float' , float                , 'Length of the reference match excluding gaps' ],
            
            [ "_meta_max_confidence_for_qry"          , 'float' , float                , 'What is the highest confidence for this QryContigID' ],
            [ "_meta_max_confidence_for_qry_chrom"    , 'float' , float                , 'Which RefContigID is the highest confidence for this QryContigID' ],
        
            [ "_meta_num_orientations"                , 'int'   , int                  , 'Number of orientations for this QryContigID' ],
            [ "_meta_num_qry_matches"                 , 'int'   , int                  , 'Number of RefContigID matches for this QryContigID' ],
            [ "_meta_qry_matches"                     , 'string', str                  , 'Which chromosomes this QryContigID matches'],
            
            [ "_meta_proportion_query_len_gapped"     , 'float' , float                , '_meta_len_qry_match_gapped / QryLen' ],
            [ "_meta_proportion_query_len_no_gap"     , 'float' , float                , '_meta_len_qry_match_no_gap / QryLen' ],
            [ "_meta_proportion_sizes_gapped"         , 'float' , float                , '_meta_len_ref_match_gapped / _meta_len_qry_match_gapped' ],
            [ "_meta_proportion_sizes_no_gap"         , 'float' , float                , '_meta_len_ref_match_no_gap / _meta_len_qry_match_no_gap' ]
        ]
}


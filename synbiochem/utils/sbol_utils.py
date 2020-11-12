'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=superfluous-parens
# pylint: disable=too-many-arguments
import re
import uuid
from xml.etree import ElementTree

from synbiochem.utils.dna_utils import DNA


_RDF_NS = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
_NS = {'ns': 'http://sbols.org/v1#',
       'rdf': _RDF_NS}


def read(filename):
    '''Parses SBOL v1.'''
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    dna_comp = root.find('ns:DnaComponent', _NS)
    dna_seq = dna_comp.find('ns:dnaSequence', _NS)
    dna_seq = dna_seq.find('ns:DnaSequence', _NS)

    params = _read_dna_comp(dna_comp)
    params.update({'seq': dna_seq.find('ns:nucleotides', _NS).text})
    dna = DNA(**params)

    for annot in dna_comp.findall('ns:annotation', _NS):
        _read_annot(dna, annot)

    return dna


def write(dna, filename=None):
    '''Writes a Dna object to SBOL v1.'''
    root = ElementTree.Element('ns2:RDF', {'xmlns': 'http://sbols.org/v1#',
                                           'xmlns:ns2': _RDF_NS})

    dna_comp = _write_dna_comp(root, dna)

    dna_seq = _write(dna_comp, 'dnaSequence')
    dna_seq = _write(dna_seq, 'DnaSequence', _get_about())
    _write(dna_seq, 'nucleotides', text=dna['seq'])

    for feature in dna['features']:
        annot = _write(dna_comp, 'annotation')
        annot = _write(annot, 'SequenceAnnotation', _get_about())
        _write(annot, 'bioStart', text=str(feature['start']))
        _write(annot, 'bioEnd', text=str(feature['end']))
        _write(annot, 'strand', text='+' if feature['forward'] else '-')
        sub_comp = _write(annot, 'subComponent')
        _write_dna_comp(sub_comp, feature)

    sbol = ElementTree.tostring(root)

    with open(filename, 'wb') as outfile:
        outfile.write(sbol)

    return sbol


def _read_dna_comp(dna_comp):
    '''Read DNAComponent node.'''
    disp_id = _find_text(dna_comp, 'ns:displayId')
    name = _find_text(dna_comp, 'ns:name')
    desc = _find_text(dna_comp, 'ns:description')
    typ_node = dna_comp.find('rdf:type', _NS)
    typ = typ_node.attrib['{' + _RDF_NS + '}resource'] \
        if typ_node is not None else None

    if not re.match('^[\\w-]+$', disp_id):
        disp_id = str(uuid.uuid4())

    return {'disp_id': disp_id, 'name': name, 'desc': desc, 'typ': typ}


def _read_annot(dna, annot):
    '''Reads annotation node.'''
    seq_annot = annot.find('ns:SequenceAnnotation', _NS)
    sub_comp = seq_annot.find('ns:subComponent', _NS)
    dna_comp = sub_comp.find('ns:DnaComponent', _NS)

    params = _read_dna_comp(dna_comp)

    start = _find_text(seq_annot, 'ns:bioStart')

    if start:
        params.update({'start': int(start)})

    end = _find_text(seq_annot, 'ns:bioEnd')

    if end:
        params.update({'end': int(end)})

    forward = _find_text(seq_annot, 'ns:strand')

    if forward:
        params.update({'forward': forward == '+'})

    # Tests due to ICE eccentricities...
    try:
        feat = DNA(**params)

        pos = (feat['start'], feat['end'], feat['forward'])

        # Prevents cases where features are duplicated in the SBOL:
        if pos not in [(feature['start'], feature['end'], feature['forward'])
                       for feature in dna['features']]:
            dna['features'].append(feat)

    except ValueError:
        # Prevents cases with no end position and no sequence:
        print('Ignoring invalid feature.')


def _find_text(parent, field):
    '''Finds text from node.'''
    node = parent.find(field, _NS)
    return None if node is None else node.text


def _write_dna_comp(parent, dna):
    '''Write DNAComponent node.'''
    dna_comp = ElementTree.SubElement(parent, 'DnaComponent',
                                      _get_about(dna['disp_id']))

    if dna['typ']:
        _write(dna_comp, 'ns2:type', {'ns2:resource': dna['typ']})

    _write(dna_comp, 'displayId', text=dna['disp_id'])

    if dna['name']:
        _write(dna_comp, 'name', text=dna['name'])

    if dna['desc']:
        _write(dna_comp, 'description', text=dna['desc'])

    return dna_comp


def _write(parent, name, params=None, text=None):
    '''Writes a node.'''
    if params is None:
        params = {}

    node = ElementTree.SubElement(parent, name, params)
    node.text = text
    return node


def _get_about(uid=None):
    '''Gets about attributes.'''
    if uid is None:
        uid = str(uuid.uuid4())

    return {'ns2:about': 'https://www.synbiochem.co.uk#' + uid}

def writeGenbank(dna, filename=None):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    #print(dna)
    # Create a sequence
    sequence_string = dna['seq']
    sequence_object = Seq(sequence_string, IUPAC.unambiguous_dna)
    # Create a record
    record = SeqRecord(sequence_object,
                       id='123456789', # random accession number
                       name='partsgenie',
                       description=dna['name'])
    # Add annotation
    for feature in dna['features']:
        print(feature['typ'])
        if feature['typ'] == "http://identifiers.org/so/SO:0000316":
             feature['typ'] ="CDS"
        elif feature['typ'] == None: print("balls")
        else:
             try:
                 feature['typ'] =feature['name']
                 try:
                    if len(feature['typ']) > 15:
                        feature['typ']  = "sequence"
                 except:
                        feature['typ']  = "sequence"
             except AssertionError:
                 feature['typ'] ="sequence"
                 print("error")
        #feature2 = SeqFeature(location=FeatureLocation(start=feature['start'], end=feature['end'] ), id=feature['name'], type=feature['typ'],  strand= 1 if feature['forward'] else -1)
        try:
             location= FeatureLocation(start=feature['start'], end=feature['end'] , strand= 1 if feature['forward'] else -1 )
        except:
             print("location")
       
        feature2 = SeqFeature(location,  id=feature['name'], type=feature['typ']
)
        try:
             feature2.qualifiers["label"] = feature['name']
        except:
             print("label")
        try:
             feature2.qualifiers["translation"]= feature['seq']
        except:
             print("transcaltion")
        try:
             feature2.qualifiers["strand"] =  strand= 1 if feature['forward'] else -1
        except:
             print("strand")
        try:
            print(feature2)
            #if feature2['typ'] ==None: feature2['typ']="blank"
            record.features.append(feature2)
        except:
            print("pain")

    output_file = open(filename, 'w')
    SeqIO.write(record, output_file, 'genbank')










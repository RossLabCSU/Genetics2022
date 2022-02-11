
from goatools import obo_parser

def main():

    protein = 'Q8N5F7'
    h = open('Hsapiens.gaf')
    output = open('TableS12_' + protein + '_GO_Annotations.tsv', 'w')
    header = header = ['Database', 'Protein ID', 'Common Gene Name', 'Qualifier', 'GO id', 'GO name', 'GO namespace', 'GO definition', 'Literature reference', 'Evidence code', '', 'Aspect (Function, Process, Component)', 'Object symbol', 'Object synonym', 'Object type', 'Tax ID', 'Date of annotation', 'Assigned by']
    output.write('\t'.join(header) + '\n')
    go_descriptions = get_go_descriptions()
    
    for line in h:
        if line.startswith('!'):
            continue
        
        items = line.rstrip().split('\t')
        if items[1] != protein:
            continue
        else:
            go_id = items[4]
            if go_id in go_descriptions:
                items = items[:5] + [go_descriptions[go_id]['name'], go_descriptions[go_id]['namespace'], go_descriptions[go_id]['definition']] + items[5:]
            else:
                items = items[:5] + ['not available']*3 + items[5:]
            output.write('\t'.join(items) + '\n')

    h.close()
    
    
def get_go_descriptions():
    
    df = {}
    h = open('go-basic.obo')
    line = ''
    while line != '[Term]':
        line = h.readline().rstrip()
        
    id = ''
    name = ''
    namespace = ''
    definition = ''
    for line in h:
        if line.rstrip() == '[Term]':
            df[id] = {'name':name, 'namespace':namespace, 'definition':definition}
            id = ''
            name = ''
            namespace = ''
            definition = ''

        if line.startswith('id: '):
            id = line.rstrip().replace('id: ', '')
        if line.startswith('name: '):
            name = line.rstrip().replace('name: ', '')
        if line.startswith('namespace: '):
            namespace = line.rstrip().replace('namespace: ', '')
        if line.startswith('def: '):
            definition = line.rstrip().replace('def: ', '')
            
    df[id] = {'name':name, 'namespace':namespace, 'definition':definition}
    
    h.close()

    return df


if __name__ == '__main__':
    main()
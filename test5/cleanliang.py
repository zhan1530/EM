import collections
import re
import json
'''

'''


testall_path = './aln.sam'
fill1_path = './fill1all'
fill2_path = './fill2all'
fill3_path = './fill3all'
fill4_path = './fill4all'
fill5_path = './fill5all'

NA = 'NA'


def stage1():
    global testall_path,fill1_path,fill2_path

    fill1_fp = open(fill1_path, 'w')
    fill2_fp = open(fill2_path, 'w')

    fill1_header = '\t'.join(['Name','V4','V6','V6m','V6t','V8','V9','L/R']) + '\n'
    fill1_fp.write(fill1_header)
    fill2_header = '\t'.join(['Name', 'V2', 'V4', 'V6', 'V6m', 'V6t', 'V8', 'V9', 'L/R']) + '\n'
    fill2_fp.write(fill2_header)

    v6_regex = re.compile(r'(\d+)(\D)')

   
    stage1_merged = {}

    last_LR = None 

    with open(testall_path, 'r') as fp:
        for row in fp:
            row = row.strip()
            if not row.startswith('M00302'):
                continue
            fields = row.split('\t')

            v6_split = v6_regex.findall(fields[5])
            item = collections.OrderedDict({
                'Name': fields[0],
                'V2': NA,
                'V4': fields[3],
                'V6': fields[5],
                'V6m': str(calc_v6m(v6_split)),
                'V6t': str(calc_v6t(v6_split)),
                'V8': fields[7],
                'V9': fields[8],
            })

            is_supply = False
            if float(fields[1]) > 2000: 
                is_supply = True
                item['V2'] = fields[1]
                item['L/R'] = last_LR
            else:
                del item['V2']
                if int(fields[8]) == 0:
                    item['L/R'] = 'L' if int(fields[3]) <= int(fields[7]) else 'R'
                else:
                    item['L/R'] = 'L' if int(fields[8]) >= 0 else 'R'

            output_row = '\t'.join(item.values()) + '\n'
            if not is_supply:
                fill1_fp.write(output_row)
            else:
                fill2_fp.write(output_row)

         
            stage1_merged.setdefault(item['Name'], {})
            sub_key = item['L/R']
            if is_supply:
                sub_key += '_S'
            stage1_merged[item['Name']][sub_key] = item

            last_LR = item['L/R']
    return stage1_merged


def stage2(stage1_output):
    global fill3_path, fill4_path

    fill3_fp = open(fill3_path, 'w')
    fill4_fp = open(fill4_path, 'w')
    fill5_fp = open(fill5_path, 'w')

    fill3_header = '\t'.join(['Name','L_V4_1','L_V6m_1','L_V6_1','L_V4_s','L_V6m_s','L_V6_s','L_V2_s','L_V6t', 'R_V4_1','R_V6m_1','R_V6_1','R_V4_s','R_V6m_s','R_V6_s', 'R_V2_s','R_V6t','V9']) + '\n'
    fill3_fp.write(fill3_header)
    fill4_header = '\t'.join(['Name','L_V4_1','L_V6m_1','L_V6_1','L_V4_2','L_V6m_2','L_V6_2','L_V2_s','L_V6t', 'R_V4_1','R_V6m_1','R_V6_1','R_V4_2','R_V6m_2','R_V6_2', 'R_V2_s','R_V6t','V9']) + '\n'
    fill4_fp.write(fill4_header)
    fill5_header = '\t'.join(['Name','L_V4_1','L_V6m_1','L_V6_1','L_V4_2','L_V6m_2','L_V6_2','L_V2_s','L_V6t', 'R_V4_1','R_V6m_1','R_V6_1','R_V4_2','R_V6m_2','R_V6_2', 'R_V2_s','R_V6t','V9']) + '\n'
    fill5_fp.write(fill5_header)

    for name, item in stage1_output.items():
        L = item.get('L')
        R = item.get('R')
        L_S = item.get('L_S')
        R_S = item.get('R_S')
        item = collections.OrderedDict({
            'Name': name,
            'L_V4_1': L['V4'] if L else NA,
            'L_V6m_1': L['V6m'] if L else NA,
            'L_V6_1': L['V6'] if L else NA,
            'L_V4_s': L_S['V4'] if L_S else NA,
            'L_V6m_s': L_S['V6m'] if L_S else NA,
            'L_V6_s': L_S['V6'] if L_S else NA,
            'L_V2_s': L_S['V2'] if L_S else NA,
            'L_V6t': L['V6t'] if L else NA,
            'R_V4_1': R['V4'] if R else NA,
            'R_V6m_1': R['V6m'] if R else NA,
            'R_V6_1': R['V6'] if R else NA,
            'R_V4_s': R_S['V4'] if R_S else NA,
            'R_V6m_s': R_S['V6m'] if R_S else NA,
            'R_V6_s': R_S['V6'] if R_S else NA,
            'R_V2_s': R_S['V2'] if R_S else NA,
            'R_V6t': R['V6t'] if R else NA,
            'V9': L['V9'] if L else NA,
        })

       
        output_row = '\t'.join(item.values()) + '\n'
        fill3_fp.write(output_row)

      
        fill4_item = {
            'Name': name
        }
        if item['L_V4_s'] == NA or int(item['L_V4_1']) < int(item['L_V4_s']):
            fill4_item['L_V4_1'] = item['L_V4_1']
            fill4_item['L_V6m_1'] = item['L_V6m_1']
            fill4_item['L_V6_1'] = item['L_V6_1']

            fill4_item['L_V4_2'] = item['L_V4_s']
            fill4_item['L_V6m_2'] = item['L_V6m_s']
            fill4_item['L_V6_2'] = item['L_V6_s']
        else:
            fill4_item['L_V4_1'] = item['L_V4_s']
            fill4_item['L_V6m_1'] = item['L_V6m_s']
            fill4_item['L_V6_1'] = item['L_V6_s']

            fill4_item['L_V4_2'] = item['L_V4_1']
            fill4_item['L_V6m_2'] = item['L_V6m_1']
            fill4_item['L_V6_2'] = item['L_V6_1']
        fill4_item['L_V2_s'] = item['L_V2_s']
        fill4_item['L_V6t'] = item['L_V6t']

        if item['R_V4_s'] == NA or int(item['R_V4_1']) < int(item['R_V4_s']):
            fill4_item['R_V4_1'] = item['R_V4_1']
            fill4_item['R_V6m_1'] = item['R_V6m_1']
            fill4_item['R_V6_1'] = item['R_V6_1']

            fill4_item['R_V4_2'] = item['R_V4_s']
            fill4_item['R_V6m_2'] = item['R_V6m_s']
            fill4_item['R_V6_2'] = item['R_V6_s']
        else:
            fill4_item['R_V4_1'] = item['R_V4_s']
            fill4_item['R_V6m_1'] = item['R_V6m_s']
            fill4_item['R_V6_1'] = item['R_V6_s']

            fill4_item['R_V4_2'] = item['R_V4_1']
            fill4_item['R_V6m_2'] = item['R_V6m_1']
            fill4_item['R_V6_2'] = item['R_V6_1']
        fill4_item['R_V2_s'] = item['R_V2_s']
        fill4_item['R_V6t'] = item['R_V6t']

        fill4_item['V9'] = item['V9']

        output_row = '\t'.join(fill4_item.values()) + '\n'
        fill4_fp.write(output_row)

        # fill5输出
        fill5_item = fill4_item
        if fill4_item['L_V4_1'] == NA or int(fill4_item['L_V4_1']) > int(fill4_item['R_V4_1']):
            fill5_item = {
                'Name': fill4_item['Name'],
                'L_V4_1': fill4_item['R_V4_1'],
                'L_V6m_1': fill4_item['R_V6m_1'],
                'L_V6_1': fill4_item['R_V6_1'],
                'L_V4_2': fill4_item['R_V4_2'],
                'L_V6m_2': fill4_item['R_V6m_2'],
                'L_V6_2': fill4_item['R_V6_2'],
                'L_V2_s': fill4_item['R_V2_s'],
                'L_V6t': fill4_item['R_V6t'],
                'R_V4_1': fill4_item['L_V4_1'],
                'R_V6m_1': fill4_item['L_V6m_1'],
                'R_V6_1': fill4_item['L_V6_1'],
                'R_V4_2': fill4_item['L_V4_2'],
                'R_V6m_2': fill4_item['L_V6m_2'],
                'R_V6_2': fill4_item['L_V6_2'],
                'R_V2_s': fill4_item['L_V2_s'],
                'R_V6t': fill4_item['L_V6t'],
                'V9': fill4_item['V9'],
            }
        output_row = '\t'.join(fill5_item.values()) + '\n'
        fill5_fp.write(output_row)


def calc_v6m(v6_split):
    s = 0
    reserved = 0
    is_begin = False

    for item in v6_split:
        if item[1] == 'M':
            reserved = 0
            is_begin = True
        if is_begin:
            s += int(item[0])
            if item[1] != 'M': 
                reserved += int(item[0])
    s -= reserved
    return s

def calc_v6t(v6_split):
    s = 0
    for item in v6_split:
        s += int(item[0])
    return s

if __name__ == '__main__':
    s1_output = stage1()
    stage2(s1_output)

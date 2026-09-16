#!PYTHONEXE

"""checkParmsCodes.py checks a parmCodes.txt file for errors

$Id: checkParmsCodes.py 7045 2019-10-07 19:56:46Z brideout $
"""

usage = 'python checkParmsCodes.py <parmCodesFile>'

import sys

# duplicates
dup_dict = {410:[411],
            414:[415],
            420:[421],
            800:[802,803],
            805:[806],
            1240:[1241,1242],
            1250:[1251,1252],
            1270:[1272],
            1280:[1282],
            1430:[1431],
            1455:[1456]}

if len(sys.argv) != 2:
    print(usage)
    sys.exit(-1)
    
parm_codes = [] # make sure all codes are unique
mnem_list = [] # make sure all mnemonics are unique

f = open(sys.argv[1])
lines = f.readlines()
f.close()
for line in lines:
    if line[0] == '#':
        continue
    if len(line.strip()) == 0:
        continue
    items = line.split(',')
    try:
        code = int(items[0])
    except:
        raise IOError('line <%s> does not start with a code' % (line))
    if len(items) != 9:
        if len(items) != 10:
            raise IOError('Illegal number of items in line <%s>' % (line))
        else:
            main_code = int(items[9])
            if main_code not in list(dup_dict.keys()):
                raise IOError('Line <%s> is not a legal duplicate parameter' % (line))
            if code not in dup_dict[main_code]:
                raise IOerror('Line <%s> is not a legal duplicate parameter' % (line))
    if code != 0:
        if code in parm_codes:
            raise IOError('duplicate code %i in line <%s>' % (code, line))
        parm_codes.append(code)
    mnem = items[3].strip().upper()
    if mnem in mnem_list:
        raise IOError('duplicate mnem %s in line <%s>' % (mnem, line))
    mnem_list.append(mnem)
    
    format = items[4]
    try:
        s = format % (1)
    except:
        raise IOError('Illegal format <%s> in line <%s>' % (format, line))
    
    try:
        width = int(items[5])
    except:
        raise IOError('Illegal width field in line %s' % (line))
    if width < 3:
        raise IOError('Width less than 3 in line <%s>' % (line))
    
    try:
        categoryId = int(items[6])
    except:
        raise IOError('Illegal categoryId field in line %s' % (line))
    if categoryId < 0:
        raise IOError('Negative categoryId in line <%s>' % (line))
    
    try:
        hasDesc = int(items[7])
    except:
        raise IOError('Illegal hasDesc field in line %s' % (line))
    if hasDesc not in (0,1):
        raise IOError('Illegal hasDesc %i in line <%s>' % (hasDesc, line))
    
    try:
        hasErrDesc = int(items[8])
    except:
        raise IOError('Illegal hasErrDesc field in line %s' % (line))
    if hasErrDesc not in (0,1):
        raise IOError('Illegal hasErrDesc %i in line <%s>' % (hasErrDesc, line))
    
    
print('No problems found with parmCodes.txt file')
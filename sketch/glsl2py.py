text = ""
with open('font.frag', 'rt') as f:
    text = f.read()
chars = []
smoothdata = []
lineardata = []
data = []
smooth = False
inblock = False

for line in text.splitlines():
    if '/*' in line:
        line = line.split()[1]
        if line in chars: 
            smooth = True
        else: chars += [ line ]
        inblock = True
    
    if inblock and 'dlinesegment' in line:
        line = line.replace(' ','').split(',')
        first = line[1].replace('vec2(', '')
        second = line[2].replace(')', '')
        third = line[3].replace('vec2(', '')
        fourth = line[4].replace(')', '')
        
        data += [ [first,second, third, fourth] ]
        
    if inblock and '*/' in line:
        if smooth:
            #print("smooth", data)
            smoothdata += [ data ]
        else:
            #print("linear", data)
            lineardata += [ data ]
        
        data = []
            
print("chars",chars)
with open("output.py", "wt") as f:
    for i in range(len(chars)):
        print(chars[i])
        f.write('    ')
        if(i!=0):
            f.write('el')
        f.write('if char == \'' + chars[i].lower() + '\':\n')
        f.write('        lines = [');
        for entry in lineardata[i][:-1]:
            print("lin>",entry)
            f.write('['+entry[0]+','+entry[1]+','+entry[2]+','+entry[3]+'],')
        f.write('['+lineardata[i][-1][0]+','+lineardata[i][-1][1]+','+lineardata[i][-1][2]+','+lineardata[i][-1][3]+']]\n')
        f.write('        smoothlines = [');
        if len(smoothdata[i]) != 0:
            for entry in smoothdata[i][:-1]:
                print("smooth>",entry)
                f.write('['+entry[0]+','+entry[1]+','+entry[2]+','+entry[3]+'],')
            f.write('['+smoothdata[i][-1][0]+','+smoothdata[i][-1][1]+','+smoothdata[i][-1][2]+','+smoothdata[i][-1][3]+']]\n')
        
    f.close()

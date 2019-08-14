lines = None
with open("sfx.frag", "rt") as f:
    lines = f.readlines()

pos_B = []
pos_t = []
pos_SPB = []
for line in lines:
    if 'const float pos_B[' in line:
        line = line.strip().split('(')[1].replace(');','').split(',')
        pos_B = [ float(entry) for entry in line ]
    if 'const float pos_SPB[' in line:
        line = line.strip().split('(')[1].replace(');','').split(',')
        pos_SPB = [ float(entry) for entry in line ]
    if 'const float pos_t[' in line:
        line = line.strip().split('(')[1].replace(');','').split(',')
        pos_t = [ float(entry) for entry in line ]

print('#version 130')
print('uniform float iTime;')
print('void scale(out float s)')
print('{')
nb_all = 0
for i in range(len(pos_B)-1):
    t_start = pos_t[i]
    t_end = pos_t[i+1]
    spb = pos_SPB[i]/4./4.
    nb = pos_B[i]
    #nbeats = mod(iTime, 60./29.);
    #iScale = nbeats-30./29.;
    #iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));
    print('    if(iTime >= ',t_start,' && iTime < ',t_end,')')
    print('    {')
    print('        s = mod(iTime-',t_start,',',spb,')-',.5*spb,';')
    print('        s = smoothstep(',-spb/12.,',0.,s)*(1.-smoothstep(0.,',spb/4.,',s));')
    print('    }')
print('}')
    

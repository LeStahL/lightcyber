#version 130
#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0.,1e-3,clamp(x,0.,1e-3)); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sin_(float a, float p) { return sin(2. * PI * mod(a,1.) + p); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq_(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float freqC1(float note){ return 32.7 * pow(2., note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return sin(.5*PI*float(n)); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

#define NTIME 4
const float pos_B[4] = float[4](0.,32.,40.,81.);
const float pos_t[4] = float[4](0.,68.571429,84.571429,159.116883);
const float pos_BPS[3] = float[3](.46666666375000004,.5,.5500000040243902);
const float pos_SPB[3] = float[3](2.14285715625,2.,1.818181804878049);
float BPS, SPB, BT;

const float Fsample = 44100.; // PRODUCTION: CHANGE THIS BACK TO 44100.
const float Tsample = 1./Fsample;

const float filterthreshold = 1e-3;

//TEXCODE

float s_atan(float a) { return 2./PI * atan(a); }
float squarey(float a, float edge) { return abs(a) < edge ? a : floor(4.*a+.5)*.25; } 

float drop_phase(float time, float t1, float f0, float f1)
{
    float t = min(time, t1);
    float phi = f0*t + .5*(f1-f0)/t1*t*t;

    if(time > t1)
    {
        phi += f1 * (time - t1);
    }
    return phi;
}

float lpnoise(float t, float fq)
{
    t *= fq;
    float tt = fract(t);
    float tn = t - tt;
    return mix(pseudorandom(floor(tn) / fq), pseudorandom(floor(tn + 1.0) / fq), smoothstep(0.0, 1.0, tt));
}

float reverb_phase(float t, float amt)
{
    float r = lpnoise(t, 100.0) + 0.2*lpnoise(t, 550.0) + 0.1*lpnoise(t, 1050.0)*exp(-5.*t);
    return amt * r;
}

float env_AHDSR(float x, float L, float A, float H, float D, float S, float R)
{
    return (x<A ? x/A : x<A+H ? 1 : x<A+H+D ? (1. - (1.-S)*(x-H-A)/D) : x<=L-R ? S : x<=L ? S*(L-x)/R : 0.);
}

float env_limit_length(float x, float length, float release)
{
    return clamp(x * 1e3, 0., 1.) * clamp(1 - (x-length)/release, 0., 1.);
}

float waveshape(float s, float amt, float A, float B, float C, float D, float E)
{
    float w;
    float m = sign(s);
    s = abs(s);

    if(s<A) w = B * smoothstep(0.,A,s);
    else if(s<C) w = C + (B-C) * smoothstep(C,A,s);
    else if(s<=D) w = s;
    else if(s<=1.)
    {
        float _s = (s-D)/(1.-D);
        w = D + (E-D) * (1.5*_s*(1.-.33*_s*_s));
    }
    else return 1.;
    
    return m*mix(s,w,amt);
}

float sinshape(float x, float amt, float parts)
{
    return (1.-amt) * x + amt * sign(x) * 0.5 * (1. - cos(parts*PI*x));
}

float comp_SAW(int N, float inv_N, float PW) {return inv_N * (1. - _sin(float(N)*PW));}
float comp_TRI(int N, float inv_N, float PW) {return N % 2 == 0 ? .1 * inv_N * _sin(float(N)*PW) : inv_N * inv_N * (1. - _sin(float(N)*PW));}
float comp_SQU(int N, float inv_N, float PW) {return inv_N * (minus1hochN(N) * _sin(.5*float(N)*PW + .25) - 1.);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}

float MADD(float t, float f, float p0, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)
{
    float ret = 0.;
    float INR = keyF==1 ? 1./CO : f/CO;
    float IRESQ = keyF==1 ? 1./RES_Q : 1./(RES_Q*f);
    
    float p = f*t;
    float float_N, inv_N, comp_mix, filter_N;
    for(int N=1; N<=NMAX; N+=NINC)
    {
        float_N = float(N);
        inv_N = 1./float_N;
        comp_mix = MIX < 0. ? (MIX+1.) * comp_TRI(N,inv_N,PW)  -     MIX  * comp_SAW(N,inv_N,PW)
                 : MIX < 1. ? (1.-MIX) * comp_TRI(N,inv_N,PW)  +     MIX  * comp_SQU(N,inv_N,PW)
                            : (MIX-1.) * comp_HAE(N,inv_N,PW)  + (2.-MIX) * comp_SQU(N,inv_N,PW);

        if(abs(comp_mix) < 1e-6) continue;

        filter_N = pow(1. + pow(float_N*INR,NDECAY),-.5) + RES * exp(-pow((float_N*f-CO)*IRESQ,2.));
                
        ret += comp_mix * filter_N * (_sin_(float_N * p, p0) + _sin_(float_N * p * (1.+DET), p0));
    }
    return s_atan(ret);
}

float QFM_FB(float PH, float FB) // my guessing of feedback coefficients, FB>0 'saw', FB<0 'sq'
{
    if(FB > 0.) return abs(FB) * .8*_sin(PH + .35*_sin(PH));
    else return abs(FB) * _sin(PH + .5*PI);
}

float QFM(float t, float f, float phase, float LV1, float LV2, float LV3, float LV4, float FR1, float FR2, float FR3, float FR4, float FB1, float FB2, float FB3, float FB4, float ALGO)
{
    int iALGO = int(ALGO);
    float PH1 = FR1 * f * t + phase;
    float PH2 = FR2 * f * t + phase;
    float PH3 = FR3 * f * t + phase;
    float PH4 = FR4 * f * t + phase;
    
    float LINK41 = 0., LINK42 = 0., LINK43 = 0., LINK32 = 0., LINK31 = 0., LINK21 = 0.; 
    if(iALGO == 1)       {LINK43 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 2)  {LINK42 = 1.; LINK32 = 1.; LINK21 = 1.;}    
    else if(iALGO == 3)  {LINK41 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 4)  {LINK42 = 1.; LINK43 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 5)  {LINK41 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 6)  {LINK43 = 1.; LINK32 = 1.;}
    else if(iALGO == 7)  {LINK43 = 1.; LINK32 = 1.; LINK31 = 1.;}
    else if(iALGO == 8)  {LINK21 = 1.; LINK43 = 1.;}
    else if(iALGO == 9)  {LINK43 = 1.; LINK42 = 1.; LINK41 = 1.;}
    else if(iALGO == 10) {LINK43 = 1.; LINK42 = 1.;}
    else if(iALGO == 11) {LINK43 = 1.;}

    float OP4 = LV4 * _sin(PH4 + QFM_FB(PH4, FB4));
    float OP3 = LV3 * _sin(PH3 + QFM_FB(PH3, FB3) + LINK43*OP4);
    float OP2 = LV2 * _sin(PH2 + QFM_FB(PH2, FB2) + LINK42*OP4 + LINK32*OP3);
    float OP1 = LV1 * _sin(PH1 + QFM_FB(PH1, FB1) + LINK41*OP4 + LINK31*OP3 + LINK32*OP2);
    
    float sum = OP1;
    if(LINK21 > 0.) sum += OP2;
    if(LINK31 + LINK32 > 0.) sum += OP3;
    if(LINK41 + LINK42 + LINK43 > 0.) sum += OP4;
    
    return s_atan(sum);
}

float protokick(float t, float f_start, float f_end, float fdecay, float hold, float decay, float drive, float detune, float rev_amount, float rev_hold, float rev_decay, float rev_drive)
{
    float phi = drop_phase(t, fdecay, f_start, f_end);
    float rev_phi = phi + reverb_phase(t, 1.);
    return clamp(drive*.5*(_sin(phi)+_sin((1.-detune)*phi)),-1.,1.) * exp(-max(t-hold, 0.)/decay)
         + rev_amount*clamp(rev_drive*.5*(_sin(rev_phi)+_sin((1.-detune)*rev_phi)),-1.,1.) * exp(-max(t-rev_hold, 0.)/rev_decay);
} 

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;
uniform int iTexSize;
uniform sampler2D iSequence;
uniform float iSequenceWidth;

// Read short value from texture at index off
float rshort(float off)
{
    float hilo = mod(off, 2.);
    off *= .5;
    vec2 ind = (vec2(mod(off, iSequenceWidth), floor(off/iSequenceWidth))+.05)/iSequenceWidth;
    vec4 block = texture(iSequence, ind);
    vec2 data = mix(block.rg, block.ba, hilo);
    return round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
float rfloat(int off)
{
    float d = rshort(float(off));
    float sign = floor(d/32768.),
        exponent = floor(d/1024.-sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

#define NTRK 5
#define NMOD 29
#define NPTN 11
#define NNOT 523
#define NDRM 10

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);} // have to put that in: "predraw" - shift every B_on earlier (and call release "prolong")
float trk_pre(int index)    {return     rfloat(index+1+4*NTRK);}
float trk_slide(int index)  {return     rfloat(index+1+5*NTRK);} // idea for future: change to individual note_slide_time
float mod_on(int index)     {return     rfloat(index+1+6*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+6*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+6*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+6*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+6*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_pan(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+3*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+4*NNOT);}
float note_slide(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+5*NNOT);}
float note_aux(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+6*NNOT);}
float drum_rel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+7*NNOT);}

vec2 mainSynth(float time)
{
    float sL = 0.;
    float sR = 0.;
    float dL = 0.;
    float dR = 0.;

    time = mod(time, 160.753247);
    
    int _it;
    for(_it = 0; _it < NTIME - 2 && pos_t[_it + 1] < time; _it++);
    BPS = pos_BPS[_it];
    SPB = pos_SPB[_it];
    BT = pos_B[_it] + (time - pos_t[_it]) * BPS;
    
    float time2 = time - .002;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, pre, f, amtL, amtR, env, aux;
    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, syn, drum;

    for(int trk = 0; trk < NTRK; trk++)
    {
        tsep0 = trk_sep(trk);
        tsep1 = trk_sep(trk + 1);

        syn = trk_syn(trk);
        rel = trk_rel(trk);
        pre = trk_pre(trk);
 
        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1) - pre); _modU++);             
        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + rel); _modL++);

        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            B = BT - mod_on(_mod);

            ptn   = mod_ptn(_mod);
            psep0 = ptn_sep(ptn);
            psep1 = ptn_sep(ptn + 1);
                         
            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1) - pre); _noteU++);
            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + rel); _noteL++);

            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                if(syn == 62)
                {
                    drum = int(note_pitch(_note));
                    rel = drum_rel(drum);
                }

                amaysynL  = 0.;
                amaysynR  = 0.;
                amaydrumL = 0.;
                amaydrumR = 0.;

                Bon   = note_on(_note) - pre;
                Boff  = note_off(_note) + rel;
                L     = Boff - Bon;
                tL    = L * SPB;
                Bprog = B - Bon;
                Bproc = Bprog / L;
                _t    = Bprog * SPB;
                _t2   = _t - .002;
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.);
                aux   = note_aux(_note);

                if(syn == 62)
                {
                    env = trk_norm(trk) * theta(Bprog) * theta(L - Bprog);
                    if(drum == 0) { sidechain = min(sidechain, 1. - vel * (clamp(1e4 * Bprog,0.,1.) - pow(Bprog/(L-rel),8.))); }
                    
                    dL += amtL * s_atan(env * amaydrumL);
                    dR += amtR * s_atan(env * amaydrumR);
                }
                else
                {
                    f = freqC1(note_pitch(_note) + mod_transp(_mod));

                    if(abs(note_slide(_note)) > 1e-3) // THIS IS SLIDEY BIZ
                    {
                        float Bslide = trk_slide(trk);
                        float fac = note_slide(_note) * log(2.)/12.;
                        if (Bprog <= Bslide)
                        {
                            float help = 1. - Bprog/Bslide;
                            f *= Bslide * (fhelp(fac) - help * fhelp(fac*help*help)) / Bprog;
                        }
                        else
                        {
                            f *= 1. + (Bslide * (fhelp(fac)-1.)) / Bprog;
                        }
                    }

                    env = theta(Bprog) * (1. - smoothstep(Boff-rel, Boff, B));
                    if(syn == 0){amaysynL = _sin(f*_t); amaysynR = _sin(f*_t2);}
                    else if(syn == 21){
                        amaysynL = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t+0.,6000.+200.*note_pitch(_note)));
                        amaysynR = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t2+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t2)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t2+0.,6000.+200.*note_pitch(_note)));
env = theta(Bprog)*pow(1.-smoothstep(Boff-rel, Boff, B),2);
                    }
                    else if(syn == 22){
                        amaysynL = env_AHDSR(Bprog,L,.2,0.,.1,1.,.9)*(1.0*(MADD((_t-SPB*0.000),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t-SPB*0.000),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*0.000),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*0.000),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*0.000),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*0.000),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)))
      +3.8e-01*(MADD((_t-SPB*6.400e-02),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t-SPB*6.400e-02),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*6.400e-02),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*6.400e-02),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*6.400e-02),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*6.400e-02),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)))
      +1.4e-01*(MADD((_t-SPB*1.280e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t-SPB*1.280e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.280e-01),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.280e-01),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.280e-01),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.280e-01),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)))
      +5.5e-02*(MADD((_t-SPB*1.920e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t-SPB*1.920e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.920e-01),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.920e-01),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.920e-01),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t-SPB*1.920e-01),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0))));
                        amaysynR = env_AHDSR(Bprog,L,.2,0.,.1,1.,.9)*(1.0*(MADD((_t2-SPB*0.000),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t2-SPB*0.000),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*0.000),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*0.000),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*0.000),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*0.000),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-0.000)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)))
      +3.8e-01*(MADD((_t2-SPB*6.400e-02),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t2-SPB*6.400e-02),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*6.400e-02),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*6.400e-02),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*6.400e-02),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*6.400e-02),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-6.400e-02)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)))
      +1.4e-01*(MADD((_t2-SPB*1.280e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t2-SPB*1.280e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.280e-01),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.280e-01),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.280e-01),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.280e-01),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.280e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)))
      +5.5e-02*(MADD((_t2-SPB*1.920e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)*.3+.5*s_atan(MADD((_t2-SPB*1.920e-01),f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.920e-01),1.01*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.920e-01),2.005*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.920e-01),4.02*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0)+MADD((_t2-SPB*1.920e-01),.49*f,0.,32,1,(-1.+(.6*Bproc)),(100.+(600.*(Bprog-1.920e-01)))+300.*(.5+(.5*_sin(.33*Bprog))),100.,5.,100.,.001,0.,0))));
                    }
                    else if(syn == 57){
                        amaysynL = ((QFM((_t-0.0*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-4.0e-04*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-8.0e-04*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-1.2e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-1.6e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-2.0e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-2.4e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-2.8e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-3.2e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t-3.6e-03*(1.+.5*_sin(1.1*_t))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.))*theta(Bprog)*exp(-3.*vel*Bprog)+.6*theta(Bprog)*exp(-3.*vel*Bprog)*sinshape(.8*MADD(_t,.5*f,0.,64,1,-.5,200.,10.,0.,3.,.013,.02*aux,0)+.2*MADD(_t,.999*f,0.,64,1,-.7,500.,10.,0.,3.,.013,.01*aux,0),.01*aux,3.));
                        amaysynR = ((QFM((_t2-0.0*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-4.0e-04*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-8.0e-04*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-1.2e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-1.6e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-2.0e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-2.4e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-2.8e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-3.2e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.)+QFM((_t2-3.6e-03*(1.+.5*_sin(1.1*_t2))),f,0.,.00787*101.,.00787*env_AHDSR(Bprog,L,.122*vel,.138,.205,.118,0.)*117.,.00787*env_AHDSR(Bprog,L,.228*vel,.249,.243,.338,0.)*41.,.00787*env_AHDSR(Bprog,L,.145*vel,.121,.183,.268,0.)*29.,.5,1.,1.001,1.,.00787*84.,.00787*31.,.00787*87.,.00787*122.,4.))*theta(Bprog)*exp(-3.*vel*Bprog)+.6*theta(Bprog)*exp(-3.*vel*Bprog)*sinshape(.8*MADD(_t2,.5*f,0.,64,1,-.5,200.,10.,0.,3.,.013,.02*aux,0)+.2*MADD(_t2,.999*f,0.,64,1,-.7,500.,10.,0.,3.,.013,.01*aux,0),.01*aux,3.));
                    }
                    else if(syn == 60){
                        amaysynL = ((QFM(floor((11000-110.*aux)*(_t-0.0*(1.+3.*_sin(.1*_t)))+.5)/(11000-110.*aux),f,0.,.00787*116.,.00787*env_AHDSR(Bprog,L,.066*vel,.037,.212,.232,0.)*9.,.00787*env_AHDSR(Bprog,L,.222*vel,.099,.224,.231,0.)*57.,.00787*env_AHDSR(Bprog,L,.206*vel,.184,.232,.243,0.)*11.,.5,1.,1.001,1.,.00787*5.,.00787*13.,.00787*98.,.00787*88.,4.)+QFM(floor((11000-110.*aux)*(_t-4.0e-03*(1.+3.*_sin(.1*_t)))+.5)/(11000-110.*aux),f,0.,.00787*116.,.00787*env_AHDSR(Bprog,L,.066*vel,.037,.212,.232,0.)*9.,.00787*env_AHDSR(Bprog,L,.222*vel,.099,.224,.231,0.)*57.,.00787*env_AHDSR(Bprog,L,.206*vel,.184,.232,.243,0.)*11.,.5,1.,1.001,1.,.00787*5.,.00787*13.,.00787*98.,.00787*88.,4.)+QFM(floor((11000-110.*aux)*(_t-8.0e-03*(1.+3.*_sin(.1*_t)))+.5)/(11000-110.*aux),f,0.,.00787*116.,.00787*env_AHDSR(Bprog,L,.066*vel,.037,.212,.232,0.)*9.,.00787*env_AHDSR(Bprog,L,.222*vel,.099,.224,.231,0.)*57.,.00787*env_AHDSR(Bprog,L,.206*vel,.184,.232,.243,0.)*11.,.5,1.,1.001,1.,.00787*5.,.00787*13.,.00787*98.,.00787*88.,4.))*env_AHDSR(Bprog,L,.099*vel,.047,.194,.289,.182));
                        amaysynR = ((QFM(floor((11000-110.*aux)*(_t2-0.0*(1.+3.*_sin(.1*_t2)))+.5)/(11000-110.*aux),f,0.,.00787*116.,.00787*env_AHDSR(Bprog,L,.066*vel,.037,.212,.232,0.)*9.,.00787*env_AHDSR(Bprog,L,.222*vel,.099,.224,.231,0.)*57.,.00787*env_AHDSR(Bprog,L,.206*vel,.184,.232,.243,0.)*11.,.5,1.,1.001,1.,.00787*5.,.00787*13.,.00787*98.,.00787*88.,4.)+QFM(floor((11000-110.*aux)*(_t2-4.0e-03*(1.+3.*_sin(.1*_t2)))+.5)/(11000-110.*aux),f,0.,.00787*116.,.00787*env_AHDSR(Bprog,L,.066*vel,.037,.212,.232,0.)*9.,.00787*env_AHDSR(Bprog,L,.222*vel,.099,.224,.231,0.)*57.,.00787*env_AHDSR(Bprog,L,.206*vel,.184,.232,.243,0.)*11.,.5,1.,1.001,1.,.00787*5.,.00787*13.,.00787*98.,.00787*88.,4.)+QFM(floor((11000-110.*aux)*(_t2-8.0e-03*(1.+3.*_sin(.1*_t2)))+.5)/(11000-110.*aux),f,0.,.00787*116.,.00787*env_AHDSR(Bprog,L,.066*vel,.037,.212,.232,0.)*9.,.00787*env_AHDSR(Bprog,L,.222*vel,.099,.224,.231,0.)*57.,.00787*env_AHDSR(Bprog,L,.206*vel,.184,.232,.243,0.)*11.,.5,1.,1.001,1.,.00787*5.,.00787*13.,.00787*98.,.00787*88.,4.))*env_AHDSR(Bprog,L,.099*vel,.047,.194,.289,.182));
                    }
                    else if(syn == 61){
                        amaysynL = (env_AHDSR(Bprog,L,0.,.2,.2,0.,0.)*env_AHDSR(Bprog,L,.001,0.,.1,1.,0.)*(1.0*env_limit_length((Bprog-0.000),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t-SPB*0.000),.2*(2.*fract(2.*f*(_t-SPB*0.000)+.4*_tri(.5*f*(_t-SPB*0.000)))-1.))+_sq_(1.04*.25*f*(_t-SPB*0.000),.2*(2.*fract(2.*f*(_t-SPB*0.000)+.4*_tri(.5*f*(_t-SPB*0.000)))-1.)))+.8*(2.*fract(2.*f*(_t-SPB*0.000)+.4*_tri(.5*f*(_t-SPB*0.000)))-1.)),1e-05,.15,.13,.3,.8,.8)+3.0e-01*env_limit_length((Bprog-2.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t-SPB*2.000e-01),.2*(2.*fract(2.*f*(_t-SPB*2.000e-01)+.4*_tri(.5*f*(_t-SPB*2.000e-01)))-1.))+_sq_(1.04*.25*f*(_t-SPB*2.000e-01),.2*(2.*fract(2.*f*(_t-SPB*2.000e-01)+.4*_tri(.5*f*(_t-SPB*2.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t-SPB*2.000e-01)+.4*_tri(.5*f*(_t-SPB*2.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)+9.0e-02*env_limit_length((Bprog-4.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t-SPB*4.000e-01),.2*(2.*fract(2.*f*(_t-SPB*4.000e-01)+.4*_tri(.5*f*(_t-SPB*4.000e-01)))-1.))+_sq_(1.04*.25*f*(_t-SPB*4.000e-01),.2*(2.*fract(2.*f*(_t-SPB*4.000e-01)+.4*_tri(.5*f*(_t-SPB*4.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t-SPB*4.000e-01)+.4*_tri(.5*f*(_t-SPB*4.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)+2.7e-02*env_limit_length((Bprog-6.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t-SPB*6.000e-01),.2*(2.*fract(2.*f*(_t-SPB*6.000e-01)+.4*_tri(.5*f*(_t-SPB*6.000e-01)))-1.))+_sq_(1.04*.25*f*(_t-SPB*6.000e-01),.2*(2.*fract(2.*f*(_t-SPB*6.000e-01)+.4*_tri(.5*f*(_t-SPB*6.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t-SPB*6.000e-01)+.4*_tri(.5*f*(_t-SPB*6.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)+8.1e-03*env_limit_length((Bprog-8.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t-SPB*8.000e-01),.2*(2.*fract(2.*f*(_t-SPB*8.000e-01)+.4*_tri(.5*f*(_t-SPB*8.000e-01)))-1.))+_sq_(1.04*.25*f*(_t-SPB*8.000e-01),.2*(2.*fract(2.*f*(_t-SPB*8.000e-01)+.4*_tri(.5*f*(_t-SPB*8.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t-SPB*8.000e-01)+.4*_tri(.5*f*(_t-SPB*8.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)));
                        amaysynR = (env_AHDSR(Bprog,L,0.,.2,.2,0.,0.)*env_AHDSR(Bprog,L,.001,0.,.1,1.,0.)*(1.0*env_limit_length((Bprog-0.000),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t2-SPB*0.000),.2*(2.*fract(2.*f*(_t2-SPB*0.000)+.4*_tri(.5*f*(_t2-SPB*0.000)))-1.))+_sq_(1.04*.25*f*(_t2-SPB*0.000),.2*(2.*fract(2.*f*(_t2-SPB*0.000)+.4*_tri(.5*f*(_t2-SPB*0.000)))-1.)))+.8*(2.*fract(2.*f*(_t2-SPB*0.000)+.4*_tri(.5*f*(_t2-SPB*0.000)))-1.)),1e-05,.15,.13,.3,.8,.8)+3.0e-01*env_limit_length((Bprog-2.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t2-SPB*2.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*2.000e-01)+.4*_tri(.5*f*(_t2-SPB*2.000e-01)))-1.))+_sq_(1.04*.25*f*(_t2-SPB*2.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*2.000e-01)+.4*_tri(.5*f*(_t2-SPB*2.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t2-SPB*2.000e-01)+.4*_tri(.5*f*(_t2-SPB*2.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)+9.0e-02*env_limit_length((Bprog-4.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t2-SPB*4.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*4.000e-01)+.4*_tri(.5*f*(_t2-SPB*4.000e-01)))-1.))+_sq_(1.04*.25*f*(_t2-SPB*4.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*4.000e-01)+.4*_tri(.5*f*(_t2-SPB*4.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t2-SPB*4.000e-01)+.4*_tri(.5*f*(_t2-SPB*4.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)+2.7e-02*env_limit_length((Bprog-6.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t2-SPB*6.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*6.000e-01)+.4*_tri(.5*f*(_t2-SPB*6.000e-01)))-1.))+_sq_(1.04*.25*f*(_t2-SPB*6.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*6.000e-01)+.4*_tri(.5*f*(_t2-SPB*6.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t2-SPB*6.000e-01)+.4*_tri(.5*f*(_t2-SPB*6.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)+8.1e-03*env_limit_length((Bprog-8.000e-01),1.*(L-rel),.07)*waveshape((s_atan(_sq_(.25*f*(_t2-SPB*8.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*8.000e-01)+.4*_tri(.5*f*(_t2-SPB*8.000e-01)))-1.))+_sq_(1.04*.25*f*(_t2-SPB*8.000e-01),.2*(2.*fract(2.*f*(_t2-SPB*8.000e-01)+.4*_tri(.5*f*(_t2-SPB*8.000e-01)))-1.)))+.8*(2.*fract(2.*f*(_t2-SPB*8.000e-01)+.4*_tri(.5*f*(_t2-SPB*8.000e-01)))-1.)),1e-05,.15,.13,.3,.8,.8)));
                    }
                    
                    sL += amtL * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynL);
                    sR += amtR * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynR);
                }
            }
        }
    }
    return vec2(s_atan(sidechain * sL + dL), s_atan(sidechain * sR + dR));
}

void main()
{
   float t = (iBlockOffset + (gl_FragCoord.x - .5) + (gl_FragCoord.y - .5)*iTexSize)/iSampleRate;
   vec2 y = mainSynth(t);
   vec2 v  = floor((0.5+0.5*y)*65535.0);
   vec2 vl = mod(v,256.0)/255.0;
   vec2 vh = floor(v/256.0)/255.0;
   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}

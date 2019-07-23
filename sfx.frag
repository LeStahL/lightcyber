#version 130
#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0.,1e-3,clamp(x,0.,1e-3)); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sin_(float a, float p) { return sin(2. * PI * mod(a,1.) + p); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq(float a) { return clip(50.*_sin(a)); }
float _psq_(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float freqC1(float note){ return 32.7 * pow(2., note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return round(sin(.5*PI*float(n))); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

const float BPM = 29.;
const float BPS = BPM/60.;
const float SPB = 60./BPM;

const float Fsample = 44100.; // PRODUCTION: CHANGE THIS BACK TO 44100.
const float Tsample = 1./Fsample;

const float stereo_delay = 2e-4; //enhance the stereo feel - this is experimental since I included the stereo functionality

const float filterthreshold = 1e-3;

//TEXCODE

float s_atan(float a) { return 2./PI * atan(a); }
float squarey(float a, float edge) { return abs(a) < edge ? a : floor(4.*a+.5)*.25; } 

float doubleslope(float t, float a, float d, float s)
{
    return smoothstep(-.00001,a,t) - (1.-s) * smoothstep(0.,d,t-a);
}

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
    float att = x/A;
    float dec = 1. - (1.-S)*(x-H-A)/D;
    float rel = (x <= L-R) ? 1. : (L-x)/R;
    return (x<A ? att : x<A+H ? 1 : x<A+H+D ? dec : x<=L-R ? S : x<=L ? (L-x)/R : 0.);
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


float comp_SAW(int N, float inv_N) {return inv_N * minus1hochN(N);}
float comp_TRI(int N, float inv_N) {return N % 2 == 0 ? 0. : inv_N * inv_N * minus1hochNminus1halbe(N);}
float comp_SQU(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (minus1hochN(N)*_sin(PW*float(N)+.25) - 1.);}

float MADD(float t, float f, float p0, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)
{
    float ret = 0.;
    float INR = keyF==1 ? 1./CO : f/CO;
    float IRESQ = keyF==1 ? 1./RES_Q : 1./(RES_Q*f);
    
    float p = f*t;
    for(int N=1; N<=NMAX; N+=NINC)
    {
        float float_N = float(N);
        float inv_N = 1./float_N;
        float comp_mix = MIX < 0. ? (MIX+1.) * comp_TRI(N,inv_N)    +  (-MIX)  * comp_SAW(N,inv_N)
                       : MIX < 1. ?   MIX    * comp_TRI(N,inv_N)    + (1.-MIX) * comp_SQU(N,inv_N,PW)
                                  : (MIX-1.) * comp_HAE(N,inv_N,PW) + (2.-MIX) * comp_SQU(N,inv_N,PW);

        float filter_N = pow(1. + pow(float_N*INR,NDECAY),-.5) + RES * exp(-pow((float_N*f-CO)*IRESQ,2.));
        
        if(abs(filter_N*comp_mix) < 1e-6) break; //or is it wise to break already?
        
        ret += comp_mix * filter_N * (_sin_(float_N * p, p0) + _sin_(float_N * p * (1.+DET), p0));
    }
    return s_atan(ret);
}

float BADD(float t, float f, float p0, float MIX, float AMP, float FPEAK, float BW, float NCUT, float DET, float PW)
{
    float ret = 0.;

    float p = f*t;
    float inv_f = 1./f;
    float inv_BW = 1./BW;
    int Nmin = int((FPEAK-NCUT*BW)*inv_f);
    int Nmax = int((FPEAK+NCUT*BW)*inv_f) + 1;
    float float_N, inv_N, comp_mix, filter_N;

    for(int N = Nmin; N < Nmax; N++)
    {
        float_N = float(N);
        inv_N = 1./float_N;
        comp_mix = MIX < 0. ? (MIX+1.) * comp_TRI(N,inv_N)    +  (-MIX)  * comp_SAW(N,inv_N)
                 : MIX < 1. ?   MIX    * comp_TRI(N,inv_N)    + (1.-MIX) * comp_SQU(N,inv_N,PW)
                 : (MIX-1.) * comp_HAE(N,inv_N,PW) + (2.-MIX) * comp_SQU(N,inv_N,PW);

        filter_N = exp(-pow((float_N*f-FPEAK)*inv_BW,2.));

        ret += comp_mix * filter_N * (_sin_(float_N * p, p0) + _sin_(float_N * p * (1.+DET), p0));
    }
    return s_atan(AMP * ret);
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



float strichr10_volume(float _BEAT)
{
    if(_BEAT<0){return 0.;}
    else if(_BEAT>=8 && _BEAT<24){return (0.+((1.-0.)/(16.-0.)*((_BEAT-8.)-0.)));}
    else if(_BEAT>=64 && _BEAT<72){return (.05+((.29-.05)/(8.-0.)*((_BEAT-64.)-0.)));}
    else if(_BEAT>=72 && _BEAT<76){return (.2+(((_BEAT-72.)-0.)));}
    else{return 1.;}
}
float strichr11_volume(float _BEAT)
{
    if(_BEAT<0){return 0.;}
    else if(_BEAT>=8 && _BEAT<12){return (.4+((1.-.4)/(4.-0.)*((_BEAT-8.)-0.)));}
    else{return 1.;}
}
float FQMimp4_volume(float _BEAT)
{
    if(_BEAT<0){return 0.;}
    else if(_BEAT>=0 && _BEAT<8){return (.3+((.65-.3)/(8.-0.)*((_BEAT-0.)-0.)));}
    else if(_BEAT>=8 && _BEAT<16){return (.65+((.4-.65)/(8.-0.)*((_BEAT-8.)-0.)));}
    else{return 1.;}
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

#define NTRK 10
#define NMOD 83
#define NPTN 19
#define NNOT 1086
#define NDRM 26

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}
float trk_slide(int index)  {return     rfloat(index+1+4*NTRK);}
float mod_on(int index)     {return     rfloat(index+1+5*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+5*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+5*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+5*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+5*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_pan(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+3*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+4*NNOT);}
float note_slide(int index) {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+5*NNOT);}
float drum_rel(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+6*NNOT);}

vec2 mainSynth(float time)
{
    float max_mod_off = 77.;
    int drum_index = 82;
    
    float sL = 0.;
    float sR = 0.;
    float dL = 0.;
    float dR = 0.;

    // mod for looping
    float BT = BPS * time;
//     float BT = mod(BPS * time, max_mod_off);
//     time = SPB * BT;
    
    float time2 = time - stereo_delay;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, f, amtL, amtR, env;
    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, syn, drum;

    for(int trk = 0; trk < NTRK; trk++)
    {
        tsep0 = trk_sep(trk);
        tsep1 = trk_sep(trk + 1);

        syn = trk_syn(trk);
        rel = trk_rel(trk);
 
        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1)); _modU++);             
        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + rel); _modL++);

        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            B = BT - mod_on(_mod);

            ptn   = mod_ptn(_mod);
            psep0 = ptn_sep(ptn);
            psep1 = ptn_sep(ptn + 1);
                         
            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1)); _noteU++);
            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + rel); _noteL++);
            //here: could introduce "monosynth" mode that sets _noteL = _noteU

            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                if(syn == drum_index)
                {
                    drum = int(note_pitch(_note));
                    rel = drum_rel(drum);
                }

                amaysynL  = 0.;
                amaysynR  = 0.;
                amaydrumL = 0.;
                amaydrumR = 0.;

                Bon   = note_on(_note);
                Boff  = note_off(_note) + rel;
                L     = Boff - Bon;
                tL    = L * SPB;
                Bprog = B - Bon;
                Bproc = Bprog / L;
                _t    = Bprog * SPB;
                _t2   = _t - stereo_delay;
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.);

                if(syn == drum_index)
                {
                    env = trk_norm(trk) * theta(Bprog) * theta(L - Bprog);
                    if(drum == 0) { sidechain = min(sidechain, 1. - vel * (clamp(1e4 * Bprog,0.,1.) - pow(Bprog/(L-rel),8.))); }
                    else if(drum == 1){
                        amaydrumL = vel*.2*fract(sin(_t*100.*.9)*50000.*.9)*doubleslope(_t,.03,.1,.1);
                        amaydrumR = vel*.2*fract(sin(_t2*100.*.9)*50000.*.9)*doubleslope(_t2,.03,.1,.1);
                    }
                    else if(drum == 2){
                        amaydrumL = vel*vel*theta(Bprog)*exp(-80.*Bprog)*lpnoise(_t + 0.,11196.)
      +vel*theta(Bprog)*exp(-70.*Bprog)*(2.*fract(12345.*_t+.7*lpnoise(_t + 0.,17896.))-1.)
      +.17*theta(Bprog)*exp(-70.*Bprog)*(2.*fract(12345.*_t+.7*lpnoise(_t + 0.,17896.))-1.);
                        amaydrumR = vel*vel*theta(Bprog)*exp(-80.*Bprog)*lpnoise(_t2 + 0.,11196.)
      +vel*theta(Bprog)*exp(-70.*Bprog)*(2.*fract(12345.*_t2+.7*lpnoise(_t2 + 0.,17896.))-1.)
      +.17*theta(Bprog)*exp(-70.*Bprog)*(2.*fract(12345.*_t2+.7*lpnoise(_t2 + 0.,17896.))-1.);
                    }
                    else if(drum == 3){
                        amaydrumL = vel*.6*(.6+(.25*_psq(4.*B)))*fract(sin(_t*100.*.3)*50000.*2.)*doubleslope(_t,0.,.09,.05)*theta(Bprog)*exp(-1.3*Bprog);
                        amaydrumR = vel*.6*(.6+(.25*_psq(4.*B)))*fract(sin(_t2*100.*.3)*50000.*2.)*doubleslope(_t2,0.,.09,.05)*theta(Bprog)*exp(-1.3*Bprog);
                    }
                    else if(drum == 5){
                        amaydrumL = vel*protokick(_t,242.,55.,.036,.088,.0666,1.42,.01,.45,.1,.15,.5)
      +.66*protokick(_t,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                        amaydrumR = vel*protokick(_t2,242.,55.,.036,.088,.0666,1.42,.01,.45,.1,.15,.5)
      +.66*protokick(_t2,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                    }
                    else if(drum == 6){
                        amaydrumL = vel*protokick(_t,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                        amaydrumR = vel*protokick(_t2,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                    }
                    else if(drum == 7){
                        amaydrumL = vel*.4*env_AHDSR(Bprog,L,.09,0.,.1,1.,0.)*(1.0*env_limit_length((Bprog-0.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-0.000))*exp(-1.45*(Bprog-0.000))))*floor(5024.*(_t-SPB*0.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*0.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*0.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*0.000)+.5)/5024.))
      +3.1e-01*env_limit_length((Bprog-5.000e-01),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-5.000e-01))*exp(-1.45*(Bprog-5.000e-01))))*floor(5024.*(_t-SPB*5.000e-01)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*5.000e-01)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*5.000e-01)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*5.000e-01)+.5)/5024.))
      +9.7e-02*env_limit_length((Bprog-1.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-1.000))*exp(-1.45*(Bprog-1.000))))*floor(5024.*(_t-SPB*1.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*1.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*1.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*1.000)+.5)/5024.))
      +3.0e-02*env_limit_length((Bprog-1.500),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-1.500))*exp(-1.45*(Bprog-1.500))))*floor(5024.*(_t-SPB*1.500)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*1.500)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*1.500)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*1.500)+.5)/5024.))
      +9.4e-03*env_limit_length((Bprog-2.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-2.000))*exp(-1.45*(Bprog-2.000))))*floor(5024.*(_t-SPB*2.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*2.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*2.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*2.000)+.5)/5024.))
      +2.9e-03*env_limit_length((Bprog-2.500),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-2.500))*exp(-1.45*(Bprog-2.500))))*floor(5024.*(_t-SPB*2.500)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*2.500)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*2.500)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*2.500)+.5)/5024.))
      +9.0e-04*env_limit_length((Bprog-3.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-3.000))*exp(-1.45*(Bprog-3.000))))*floor(5024.*(_t-SPB*3.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t-SPB*3.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t-SPB*3.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t-SPB*3.000)+.5)/5024.)));
                        amaydrumR = vel*.4*env_AHDSR(Bprog,L,.09,0.,.1,1.,0.)*(1.0*env_limit_length((Bprog-0.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-0.000))*exp(-1.45*(Bprog-0.000))))*floor(5024.*(_t2-SPB*0.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*0.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*0.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*0.000)+.5)/5024.))
      +3.1e-01*env_limit_length((Bprog-5.000e-01),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-5.000e-01))*exp(-1.45*(Bprog-5.000e-01))))*floor(5024.*(_t2-SPB*5.000e-01)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*5.000e-01)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*5.000e-01)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*5.000e-01)+.5)/5024.))
      +9.7e-02*env_limit_length((Bprog-1.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-1.000))*exp(-1.45*(Bprog-1.000))))*floor(5024.*(_t2-SPB*1.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*1.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*1.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*1.000)+.5)/5024.))
      +3.0e-02*env_limit_length((Bprog-1.500),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-1.500))*exp(-1.45*(Bprog-1.500))))*floor(5024.*(_t2-SPB*1.500)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*1.500)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*1.500)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*1.500)+.5)/5024.))
      +9.4e-03*env_limit_length((Bprog-2.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-2.000))*exp(-1.45*(Bprog-2.000))))*floor(5024.*(_t2-SPB*2.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*2.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*2.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*2.000)+.5)/5024.))
      +2.9e-03*env_limit_length((Bprog-2.500),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-2.500))*exp(-1.45*(Bprog-2.500))))*floor(5024.*(_t2-SPB*2.500)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*2.500)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*2.500)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*2.500)+.5)/5024.))
      +9.0e-04*env_limit_length((Bprog-3.000),.34*(L-rel),.29)*_sin_((400.+(1111.*theta((Bprog-3.000))*exp(-1.45*(Bprog-3.000))))*floor(5024.*(_t2-SPB*3.000)+.5)/5024.,(2.*fract(30.*floor(5024.*(_t2-SPB*3.000)+.5)/5024.+.4)-1.)+clip((1.+1.)*lpnoise(floor(5024.*(_t2-SPB*3.000)+.5)/5024. + 0.,666.))+_sin(1666.*floor(5024.*(_t2-SPB*3.000)+.5)/5024.)));
                    }
                    else if(drum == 16){
                        amaydrumL = vel*11.38*lpnoise(_t, 1.+6.77*_t)*(smoothstep(0.,1e-3,_t)-smoothstep(0.,.023,_t-.018));
                        amaydrumR = vel*11.38*lpnoise(_t2, 1.+6.77*_t2)*(smoothstep(0.,1e-3,_t2)-smoothstep(0.,.023,_t2-.018));
                    }
                    else if(drum == 17){
                        amaydrumL = vel*(.3+(.7*theta(Bprog)*exp(-10.*Bprog)))*sinshape(((clamp(.72*_tri(drop_phase(_t,.07,227.,107.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t-.09))+.69*clamp(.94*_tri(drop_phase(_t,.07,227.,107.)+.69*lpnoise(_t,8709.)),-1.,1.)*exp(-13.89*_t)+.06*lpnoise(_t,18745.)*(1.-smoothstep(0.,.67,_t-.17))+.28*lpnoise(_t,7506.)*exp(-_t*12.44)+.82*lpnoise(_t,2600.)*exp(-_t*17.54))*smoothstep(0.,.004,_t)),.6,9.);
                        amaydrumR = vel*(.3+(.7*theta(Bprog)*exp(-10.*Bprog)))*sinshape(((clamp(.72*_tri(drop_phase(_t2,.07,227.,107.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t2-.09))+.69*clamp(.94*_tri(drop_phase(_t2,.07,227.,107.)+.69*lpnoise(_t2,8709.)),-1.,1.)*exp(-13.89*_t2)+.06*lpnoise(_t2,18745.)*(1.-smoothstep(0.,.67,_t2-.17))+.28*lpnoise(_t2,7506.)*exp(-_t2*12.44)+.82*lpnoise(_t2,2600.)*exp(-_t2*17.54))*smoothstep(0.,.004,_t2)),.6,9.);
                    }
                    else if(drum == 20){
                        amaydrumL = vel*((clamp(1.09*_tri(drop_phase(_t,.08,249.,77.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t-.04))+.97*clamp(.99*_tri(drop_phase(_t,.08,249.,77.)+.97*lpnoise(_t,9855.)),-1.,1.)*exp(-21.22*_t)+.03*lpnoise(_t,10655.)*(1.-smoothstep(0.,.58,_t-.81))+.71*lpnoise(_t,7520.)*exp(-_t*16.22)+.57*lpnoise(_t,4386.)*exp(-_t*29.48))*smoothstep(0.,.005,_t));
                        amaydrumR = vel*((clamp(1.09*_tri(drop_phase(_t2,.08,249.,77.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t2-.04))+.97*clamp(.99*_tri(drop_phase(_t2,.08,249.,77.)+.97*lpnoise(_t2,9855.)),-1.,1.)*exp(-21.22*_t2)+.03*lpnoise(_t2,10655.)*(1.-smoothstep(0.,.58,_t2-.81))+.71*lpnoise(_t2,7520.)*exp(-_t2*16.22)+.57*lpnoise(_t2,4386.)*exp(-_t2*29.48))*smoothstep(0.,.005,_t2));
                    }
                    else if(drum == 21){
                        amaydrumL = vel*.57*((clamp(1.54*_tri(drop_phase(_t,.06,254.,67.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t-.25))+.48*clamp(.83*_tri(drop_phase(_t,.06,254.,67.)+.48*lpnoise(_t,1022.)),-1.,1.)*exp(-18.09*_t)+.39*lpnoise(_t,14859.)*(1.-smoothstep(0.,.51,_t-.47))+.25*lpnoise(_t,2739.)*exp(-_t*22.61)+.34*lpnoise(_t,5121.)*exp(-_t*23.75))*smoothstep(0.,.003,_t));
                        amaydrumR = vel*.57*((clamp(1.54*_tri(drop_phase(_t2,.06,254.,67.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t2-.25))+.48*clamp(.83*_tri(drop_phase(_t2,.06,254.,67.)+.48*lpnoise(_t2,1022.)),-1.,1.)*exp(-18.09*_t2)+.39*lpnoise(_t2,14859.)*(1.-smoothstep(0.,.51,_t2-.47))+.25*lpnoise(_t2,2739.)*exp(-_t2*22.61)+.34*lpnoise(_t2,5121.)*exp(-_t2*23.75))*smoothstep(0.,.003,_t2));
                    }
                    else if(drum == 22){
                        amaydrumL = vel*.42*((clamp(2.24*_tri(drop_phase(_t,.08,217.,64.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t-.29))+.7*clamp(.84*_tri(drop_phase(_t,.08,217.,64.)+.7*lpnoise(_t,1936.)),-1.,1.)*exp(-23.11*_t)+.08*lpnoise(_t,5166.)*(1.-smoothstep(0.,.16,_t-.55))+.77*lpnoise(_t,6784.)*exp(-_t*29.89)+.53*lpnoise(_t,4404.)*exp(-_t*24.64))*smoothstep(0.,.002,_t));
                        amaydrumR = vel*.42*((clamp(2.24*_tri(drop_phase(_t2,.08,217.,64.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t2-.29))+.7*clamp(.84*_tri(drop_phase(_t2,.08,217.,64.)+.7*lpnoise(_t2,1936.)),-1.,1.)*exp(-23.11*_t2)+.08*lpnoise(_t2,5166.)*(1.-smoothstep(0.,.16,_t2-.55))+.77*lpnoise(_t2,6784.)*exp(-_t2*29.89)+.53*lpnoise(_t2,4404.)*exp(-_t2*24.64))*smoothstep(0.,.002,_t2));
                    }
                    else if(drum == 23){
                        amaydrumL = vel*(.57*(0.*lpnoise(_t,981.)+0.*lpnoise(_t,950.)+0.*lpnoise(_t,2937.))*(smoothstep(0.,0.,_t)-smoothstep(0.,.13,_t-.76))+_sin(drop_phase(_t,.04,607.,288.))*exp(-_t*4.5)*.46+_sin(drop_phase(_t*1076.,.04,607.,288.))*exp(-_t*9.6)*.49);
                        amaydrumR = vel*(.57*(0.*lpnoise(_t2,981.)+0.*lpnoise(_t2,950.)+0.*lpnoise(_t2,2937.))*(smoothstep(0.,0.,_t2)-smoothstep(0.,.13,_t2-.76))+_sin(drop_phase(_t2,.04,607.,288.))*exp(-_t2*4.5)*.46+_sin(drop_phase(_t2*1076.,.04,607.,288.))*exp(-_t2*9.6)*.49);
                    }
                    else if(drum == 26){
                        amaydrumL = vel*(.25*(0.*lpnoise(_t,260.)+0.*lpnoise(_t,868.)+0.*lpnoise(_t,454.))*(smoothstep(0.,.01,_t)-smoothstep(0.,.91,_t-.91))+_sin(drop_phase(_t,.02,559.,252.))*exp(-_t*5.5)*.4+_sin(drop_phase(_t*1444.,.02,559.,252.))*exp(-_t*6.9)*.3);
                        amaydrumR = vel*(.25*(0.*lpnoise(_t2,260.)+0.*lpnoise(_t2,868.)+0.*lpnoise(_t2,454.))*(smoothstep(0.,.01,_t2)-smoothstep(0.,.91,_t2-.91))+_sin(drop_phase(_t2,.02,559.,252.))*exp(-_t2*5.5)*.4+_sin(drop_phase(_t2*1444.,.02,559.,252.))*exp(-_t2*6.9)*.3);
                    }
                    
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
                    else if(syn == 13){
                        amaysynL = ((vel*(QFM((_t-0.0*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t-4.0e-03*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t-8.0e-03*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.))*env_AHDSR(_t,tL,.082,.001,.062,.521,.153)));
                        amaysynR = ((vel*(QFM((_t-0.0*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t-4.0e-03*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t-8.0e-03*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.))*env_AHDSR(_t,tL,.082,.001,.062,.521,.153)));
                    }
                    else if(syn == 72){
                        amaysynL = (.7*vel*env_AHDSR(_t,tL,.406,0.,.1,1.,.031)*waveshape(MADD(_t,f,0.,256,1,-.523,(2598.+(810.*clip((1.+.744)*_sin(.235*B)))),9.903,1.872,6.88,.004,.68,0),.68,.1,.5,.2,.5,.6)+.6*vel*env_AHDSR(_t,tL,.406,0.,.1,1.,.231)*MADD(_t,.501*f,0.,64,1,-.523,(2598.+(810.*clip((1.+.744)*_sin(.235*B)))),9.903,1.872,6.88,.004,.2*.68,0));
                        amaysynR = (.7*vel*env_AHDSR(_t2,tL,.406,0.,.1,1.,.031)*waveshape(MADD(_t2,f,0.,256,1,-.523,(2598.+(810.*clip((1.+.744)*_sin(.235*B)))),9.903,1.872,6.88,.004,.68,0),.68,.1,.5,.2,.5,.6)+.6*vel*env_AHDSR(_t2,tL,.406,0.,.1,1.,.231)*MADD(_t2,.501*f,0.,64,1,-.523,(2598.+(810.*clip((1.+.744)*_sin(.235*B)))),9.903,1.872,6.88,.004,.2*.68,0));
                    }
                    else if(syn == 75){
                        amaysynL = strichr10_volume(BT)*smoothstep(0.,1.19701*(.6-vel+1e-3),Bprog)*(env_AHDSR(_t,tL,.24,0.,.1,1.,.262)*waveshape(MADD(_t,f,0.,256,1,-.286,(3232.+(1975.*clip((1.+.775)*_sin(.485*B)))),11.902,.763,3.63,.017,-.13,0),-.13,.1,.5,.2,.5,.6)+.5*env_AHDSR(_t,tL,.24,0.,.1,1.,.262)*MADD(_t,.501*f,0.,64,1,-.286,(3232.+(1975.*clip((1.+.775)*_sin(.485*B)))),11.902,.763,3.63,.017,.2*-.13,0));
                        amaysynR = strichr10_volume(BT)*smoothstep(0.,1.19701*(.6-vel+1e-3),Bprog)*(env_AHDSR(_t2,tL,.24,0.,.1,1.,.262)*waveshape(MADD(_t2,f,0.,256,1,-.286,(3232.+(1975.*clip((1.+.775)*_sin(.485*B)))),11.902,.763,3.63,.017,-.13,0),-.13,.1,.5,.2,.5,.6)+.5*env_AHDSR(_t2,tL,.24,0.,.1,1.,.262)*MADD(_t2,.501*f,0.,64,1,-.286,(3232.+(1975.*clip((1.+.775)*_sin(.485*B)))),11.902,.763,3.63,.017,.2*-.13,0));
                    }
                    else if(syn == 76){
                        amaysynL = strichr11_volume(BT)*(vel*env_AHDSR(_t,tL,.249,0.,.1,1.,.37)*waveshape(MADD(_t,f,0.,256,1,-.383,(7365.+(1083.*clip((1.+.916)*_sin(.525*B)))),18.724,1.488,1.24,.022,-.27,0),-.27,.1,.5,.2,.5,.6)+.5*vel*env_AHDSR(_t,tL,.249,0.,.1,1.,.37)*MADD(_t,.501*f,0.,64,1,-.383,(7365.+(1083.*clip((1.+.916)*_sin(.525*B)))),18.724,1.488,1.24,.022,.2*-.27,0));
                        amaysynR = strichr11_volume(BT)*(vel*env_AHDSR(_t2,tL,.249,0.,.1,1.,.37)*waveshape(MADD(_t2,f,0.,256,1,-.383,(7365.+(1083.*clip((1.+.916)*_sin(.525*B)))),18.724,1.488,1.24,.022,-.27,0),-.27,.1,.5,.2,.5,.6)+.5*vel*env_AHDSR(_t2,tL,.249,0.,.1,1.,.37)*MADD(_t2,.501*f,0.,64,1,-.383,(7365.+(1083.*clip((1.+.916)*_sin(.525*B)))),18.724,1.488,1.24,.022,.2*-.27,0));
                    }
                    else if(syn == 78){
                        amaysynL = .3*clip((1.+.7)*_sq_(.5*f*_t,.9+-.2*(1.+(3.*theta(Bprog)*exp(-150.*Bprog)))))
      +2.*vel*BADD(_t,f*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),0.,1.4,3.,800.*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),100.,4.,.007,0.)
      +.3*clip((1.+.2)*_sq_(.499*f*_t,-.5+.3*(1.+(3.*theta(Bprog)*exp(-150.*Bprog)))))
      +2.*vel*clip((1.+.2)*_sq_(.499*f*_t,-.5+.3*(1.+(3.*theta(Bprog)*exp(-150.*Bprog)))))
      +BADD(_t,f*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),0.,1.4,3.,800.*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),100.,4.,.007,0.);
                        amaysynR = .3*clip((1.+.7)*_sq_(.5*f*_t2,.9+-.2*(1.+(3.*theta(Bprog)*exp(-150.*Bprog)))))
      +2.*vel*BADD(_t2,f*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),0.,1.4,3.,800.*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),100.,4.,.007,0.)
      +.3*clip((1.+.2)*_sq_(.499*f*_t2,-.5+.3*(1.+(3.*theta(Bprog)*exp(-150.*Bprog)))))
      +2.*vel*clip((1.+.2)*_sq_(.499*f*_t2,-.5+.3*(1.+(3.*theta(Bprog)*exp(-150.*Bprog)))))
      +BADD(_t2,f*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),0.,1.4,3.,800.*(1.+(3.*theta(Bprog)*exp(-150.*Bprog))),100.,4.,.007,0.);
                    }
                    else if(syn == 79){
                        amaysynL = .5*(1.0*env_limit_length((Bprog-0.000),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*0.000),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*0.000),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*0.000))))
      +5.3e-01*env_limit_length((Bprog-1.330e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*1.330e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*1.330e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*1.330e-01))))
      +2.8e-01*env_limit_length((Bprog-2.660e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*2.660e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*2.660e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*2.660e-01))))
      +1.5e-01*env_limit_length((Bprog-3.990e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*3.990e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*3.990e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*3.990e-01))))
      +7.9e-02*env_limit_length((Bprog-5.320e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*5.320e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*5.320e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*5.320e-01)))))
      +.5*(1.0*env_limit_length((Bprog-0.000),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*0.000),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*0.000),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*0.000))))
      +4.3e-01*env_limit_length((Bprog-8.660e-02),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*8.660e-02),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*8.660e-02),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*8.660e-02))))
      +1.9e-01*env_limit_length((Bprog-1.732e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*1.732e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*1.732e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*1.732e-01))))
      +8.2e-02*env_limit_length((Bprog-2.598e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*2.598e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*2.598e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*2.598e-01))))
      +3.6e-02*env_limit_length((Bprog-3.464e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*3.464e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*3.464e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*3.464e-01))))
      +1.6e-02*env_limit_length((Bprog-4.330e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t-SPB*4.330e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t-SPB*4.330e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t-SPB*4.330e-01)))));
                        amaysynR = .5*(1.0*env_limit_length((Bprog-0.000),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*0.000),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*0.000),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*0.000))))
      +5.3e-01*env_limit_length((Bprog-1.330e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*1.330e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*1.330e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*1.330e-01))))
      +2.8e-01*env_limit_length((Bprog-2.660e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*2.660e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*2.660e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*2.660e-01))))
      +1.5e-01*env_limit_length((Bprog-3.990e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*3.990e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*3.990e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*3.990e-01))))
      +7.9e-02*env_limit_length((Bprog-5.320e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*5.320e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*5.320e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*5.320e-01)))))
      +.5*(1.0*env_limit_length((Bprog-0.000),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*0.000),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*0.000),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*0.000))))
      +4.3e-01*env_limit_length((Bprog-8.660e-02),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*8.660e-02),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*8.660e-02),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*8.660e-02))))
      +1.9e-01*env_limit_length((Bprog-1.732e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*1.732e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*1.732e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*1.732e-01))))
      +8.2e-02*env_limit_length((Bprog-2.598e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*2.598e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*2.598e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*2.598e-01))))
      +3.6e-02*env_limit_length((Bprog-3.464e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*3.464e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*3.464e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*3.464e-01))))
      +1.6e-02*env_limit_length((Bprog-4.330e-01),.3*(L-rel),.066)*(vel*waveshape(QFM((_t2-SPB*4.330e-01),f,0.,.00787*54.,.00787*82.,.00787*86.,.00787*64.,.5,1.,1.001,1.,.00787*95.,.00787*83.,.00787*18.,.00787*83.,10.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR((_t2-SPB*4.330e-01),tL,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*(_t2-SPB*4.330e-01)))));
                    }
                    else if(syn == 80){
                        amaysynL = FQMimp4_volume(BT)*(vel*waveshape(QFM(_t,f,0.,.00787*15.,.00787*122.,.00787*39.,.00787*53.,.5,1.,1.001,1.,.00787*81.,.00787*81.,.00787*57.,.00787*87.,9.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR(Bprog,L,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*_t)));
                        amaysynR = FQMimp4_volume(BT)*(vel*waveshape(QFM(_t2,f,0.,.00787*15.,.00787*122.,.00787*39.,.00787*53.,.5,1.,1.001,1.,.00787*81.,.00787*81.,.00787*57.,.00787*87.,9.),(.5+(.5*_sin(.4*B))),.05,.46,.3,.7,.8)*env_AHDSR(Bprog,L,.2,0.,.01,1.,.05)+vel*clip((1.+2.)*_sin(f*_t2)));
                    }
                    else if(syn == 81){
                        amaysynL = (env_AHDSR(Bprog,L,.014,0.,.01,1.,.006)*sinshape(MADD(_t,.5*f,0.,256,2,1+.9*(.55+(.4*clip((1.+1.)*_sin(4.*B)))),(1344.+(712.*_sin_(2.*B,.4))),15.,7.13,7.25,.015,.4*(.55+(.4*clip((1.+1.)*_sin(4.*B)))),0),1.,3.)+.8*clip((1.+.5)*_sin(.499*f*_t))+env_AHDSR(Bprog,L,.014,0.,.01,1.,.006)*.4*_sq_(1.01*f*_t,.95));
                        amaysynR = (env_AHDSR(Bprog,L,.014,0.,.01,1.,.006)*sinshape(MADD(_t2,.5*f,0.,256,2,1+.9*(.55+(.4*clip((1.+1.)*_sin(4.*B)))),(1344.+(712.*_sin_(2.*B,.4))),15.,7.13,7.25,.015,.4*(.55+(.4*clip((1.+1.)*_sin(4.*B)))),0),1.,3.)+.8*clip((1.+.5)*_sin(.499*f*_t2))+env_AHDSR(Bprog,L,.014,0.,.01,1.,.006)*.4*_sq_(1.01*f*_t2,.95));
                    }
                    
                    sL += amtL * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynL);
                    sR += amtR * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynR);
                }
            }
        }
    }
    return vec2(s_atan(sidechain * sL + dL), s_atan(sidechain * sR + dR));
}

vec2 mainSound(float t)
{
    return mainSynth(t);
}

void main()
{
   float t = (iBlockOffset + (gl_FragCoord.x - .5) + (gl_FragCoord.y - .5)*iTexSize)/iSampleRate;
   vec2 y = mainSound( t );
   vec2 v  = floor((0.5+0.5*y)*65535.0);
   vec2 vl = mod(v,256.0)/255.0;
   vec2 vh = floor(v/256.0)/255.0;
   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}

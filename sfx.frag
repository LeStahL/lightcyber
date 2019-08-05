#version 130
#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0.,1e-3,clamp(x,0.,1e-3)); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq_(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float freqC1(float note){ return 32.7 * pow(2., note/12.); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

#define NTIME 10
const float pos_B[10] = float[10](0.,4.,5.,6.,7.,8.,11.,11.5,12.,16.);
const float pos_t[10] = float[10](0.,8.,9.875,11.639706,13.306373,14.844834,19.130548,20.012901,21.166747,34.500081);
const float pos_BPS[9] = float[9](.5,.5333333333333333,.5666666288888913,.5999998800000239,.6500002275000797,.7000000466666697,.5666666288888924,.43333339111111824,.29999998500000075);
const float pos_SPB[9] = float[9](2.,1.875,1.7647060000000003,1.6666670000000003,1.5384609999999999,1.4285713333333334,1.7647059999999968,2.307692000000003,3.3333335);
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

#define NTRK 2
#define NMOD 20
#define NPTN 2
#define NNOT 27
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

    time = mod(time, 35.833414);
    
    int _it;
    for(_it = 0; _it < NTIME - 2 && pos_t[_it + 1] < time; _it++);
    BPS = pos_BPS[_it];
    SPB = pos_SPB[_it];
    BT = pos_B[_it] + (time - pos_t[_it]) * BPS;
    
    float time2 = time - .0002;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, pre, f, amtL, amtR, env, slide, aux;
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
                if(syn == 35)
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
                _t2   = _t - .0002;
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.);
                slide = note_slide(_note);
                aux   = note_aux(_note);

                if(syn == 35)
                {
                    env = trk_norm(trk) * theta(Bprog) * theta(L - Bprog);
                    if(drum == 0) { sidechain = min(sidechain, 1. - vel * (clamp(1e4 * Bprog,0.,1.) - pow(Bprog/(L-rel),8.))); }
                    else if(drum == 8){
                        amaydrumL = vel*(vel*(clamp(1.32*_tri(drop_phase(_t,.06,308.,80.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t-.15))+.82*clamp(.49*_tri(drop_phase(_t,.06,308.,80.)+.82*lpnoise(_t,4595.)),-1.,1.)*exp(-1.97*_t)+.09*lpnoise(_t,4032.)*(1.-smoothstep(0.,.97,_t-.79))+.1*lpnoise(_t,1111.)*exp(-_t*12.69)+.6*lpnoise(_t,7795.)*exp(-_t*1.08))*smoothstep(0.,.003,_t));
                        amaydrumR = vel*(vel*(clamp(1.32*_tri(drop_phase(_t2,.06,308.,80.)),-1.,1.)*(1.-smoothstep(-1e-3,0.,_t2-.15))+.82*clamp(.49*_tri(drop_phase(_t2,.06,308.,80.)+.82*lpnoise(_t2,4595.)),-1.,1.)*exp(-1.97*_t2)+.09*lpnoise(_t2,4032.)*(1.-smoothstep(0.,.97,_t2-.79))+.1*lpnoise(_t2,1111.)*exp(-_t2*12.69)+.6*lpnoise(_t2,7795.)*exp(-_t2*1.08))*smoothstep(0.,.003,_t2));
                    }
                    
                    dL += amtL * s_atan(env * amaydrumL);
                    dR += amtR * s_atan(env * amaydrumR);
                }
                else
                {
                    f = freqC1(note_pitch(_note) + mod_transp(_mod));

                    if(abs(slide) > 1e-3) // THIS IS SLIDEY BIZ
                    {
                        float Bslide = trk_slide(trk);
                        float fac = slide * log(2.)/12.;
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
                    else if(syn == 32){
                        amaysynL = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t+0.,6000.+200.*note_pitch(_note)));
                        amaysynR = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t2+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t2)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t2+0.,6000.+200.*note_pitch(_note)));
env = theta(Bprog)*pow(1.-smoothstep(Boff-rel, Boff, B),2);
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

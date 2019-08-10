// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float size = .005;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void dvoronoi(in vec2 x, out float d, out vec2 z)
{
    vec2 y = floor(x);
       float ret = 1.;
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    d = ret;
    z = pf;
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void scene(in vec3 x, out vec2 sdf)
{
    /*
    vec3 y = vec3(mod(x.xy,size)-.5*size, x.z);
    vec2 yi = x.xy-y.xy;
	*/
    
    float v;
    vec2 vi;
    dvoronoi(x.xy/size, v, vi);
    vec3 y = vec3(x.xy-vi*size, x.z);
    vec2 yi = vi*size;
    
    float n;
    lfnoise(4.*(yi-.5*iTime), n);
    lfnoise(12.*vec2(n,1.)*yi-1.e-2*(.8+.2*n)*iTime*c.xy, n);
    sdf = vec2(length(y-.05*n*c.yyx)-.5*size, 1.);
}

void normal(in vec3 x, out vec3 n, in float dx)
{
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = (fragCoord)/iResolution.xy*vec2(a,1.),
        s;
    
    vec3 col = c.yyy, 
        o = c.yyx+.5*vec3(cos(iTime), sin(iTime),0.),
        r = c.xyy,
        u = c.yxy,
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 300,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.05-.5*size)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        
        if(x.z<-.05-.5*size)
        {
            col = c.yyy;
            i = N;
            break;
        }
        d += min(s.x,1.e-3);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n, 5.e-4);
        vec3 l = normalize(x+.5*n);
       
		if(s.y == 1.)
        {
            float v;
            vec2 vi;
            dvoronoi(x.xy/size, v, vi);
            vec3 y = vec3(x.xy-vi*size, x.z);
            vec2 yi = vi*size;
            
            col = texture(iChannel0, yi/vec2(a,1.)).rgb;
            
            col = .4*col
                + .9*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
        }
    }
    
    //col = texture(iChannel0, uv/vec2(a,1.)).rgb;
    
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

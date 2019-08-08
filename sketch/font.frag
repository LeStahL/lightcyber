// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float iScale, nbeats;


void dbox(in vec2 x, in vec2 b, out float d)
{
    vec2 da = abs(x)-b;
    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    nbeats = mod(iTime, 60./29.);
    iScale = nbeats-30./29.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    uv *= 4.;
    
    vec3 col = vec3(0.20,0.01,0.14);
    float d, da, db;
    dbox(uv, c.xx, d);
    stroke(d, .001, d);
    col = mix(col, c.xxy, sm(d));
    
    d = 1.;
    
    /* A
    dlinesegment(uv, vec2(-.65,-.75), vec2(-.65,.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.65,0.), vec2(.65,.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.65,.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.65,0.), vec2(.65, 0.), da);
    d = min(d, da);
    */
    
    /* B
    dlinesegment(uv, vec2(-.75,.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.75, 0.), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.65,0.), vec2(.65,.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.75,-.75), da);
    d = min(d, da);
    */ 
    
    /* C
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,.35), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.35), vec2(.75, -.75), da);
    d = min(d, da);
	*/
    
    /* D
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
	*/
    
    /* E
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(.65, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.65, 0.), da);
    d = min(d, da);
    */
    
    /* F
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.65, 0.), da);
    d = min(d, da);
	*/
    
    /* G
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,.5), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.25, 0.), da);
    d = min(d, da);
	*/
    
    /* H
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.65, 0.), da);
    d = min(d, da);
    */
    
    /* I
    dlinesegment(uv, vec2(0.,-.75), vec2(0., .75), da);
    d = min(d, da);
	*/
    
    /* J
    dlinesegment(uv, vec2(.65,-.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, -.35), da);
    d = min(d, da);
    */
    
    /* K
    dlinesegment(uv, vec2(-.75,0.), vec2(.75, 0.), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.65,0.), vec2(.65,.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.75,-.75), da);
    d = min(d, da);
    */
    
    /* L
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(.65, -.75), da);
    d = min(d, da);
    */
    
    /* M
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(0.,-.75), vec2(0., .75), da);
    d = min(d, da);
	*/
    
    /* N
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
	*/
    
    /* O
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
    */
    
    /* P
    dlinesegment(uv, vec2(-.75,.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.65, 0.), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.65,0.), vec2(.65,.75), da);
    d = min(d, da);
	*/
    
    /* Q
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
    */
    
    /* R
    dlinesegment(uv, vec2(-.75,.75), vec2(.65, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.75, 0.), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.65,0.), vec2(.65,.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.75,-.75), da);
    d = min(d, da);
    */
    
    /* S
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,.5), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.5), vec2(-.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.75, 0.), da);
    d = min(d, da);
	*/
    
    /* T
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(0.,-.75), vec2(0., .75), da);
    d = min(d, da);
    */
    
    /* U
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
    */
    
    /* V
    dlinesegment(uv, vec2(-.65,.75), vec2(-.65,-.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.65,0.), vec2(.65,-.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.65,-.75), vec2(.65,-.75), da);
    d = min(d, da);
    */
    
    /* W
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(0.,-.75), vec2(0., .75), da);
    d = min(d, da);
    */
    
    /* X
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75, -.1), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.75), vec2(.75, -.1), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(-.75, .1), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,.75), vec2(.75, .1), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.65,0.), vec2(.65, 0.), da);
    d = min(d, da);
    */
    
    /* Y
    dlinesegment(uv, vec2(0.,-.75), vec2(0., 0.), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(.75, 0.), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,0.), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,0.), vec2(.75, .75), da);
    d = min(d, da);
    */
    
    /* Z
    dlinesegment(uv, vec2(-.75,-.75), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.75), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,.25), vec2(.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,.5), vec2(-.75, .75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(.75,-.5), vec2(.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.25), vec2(-.75, -.75), da);
    d = min(d, da);
    dlinesegment(uv, vec2(-.75,-.25), vec2(.75, .25), da);
    d = min(d, da);
    */
    
    
    
    stroke(d,.1,d);
    
    da = 1.;
    
    /* A
    dlinesegment(uv, vec2(-.75,-.75), vec2(-.75,0.), db);
    da = min(da, db);
    dlinesegment(uv, vec2(.75,0.), vec2(.75,-.75), db);
    da = min(da, db);
    */
    
    /* B
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* C
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* D
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
	*/
    
    /* E
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* F
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */
    
    /* G
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* H
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */
    
    /* I
    dlinesegment(uv, vec2(.1,-.75), vec2(.1,0.), db);
    da = min(da, db);
    */
    
    /* J
    dlinesegment(uv, vec2(.75,-.65), vec2(.75,0.), db);
    da = min(da, db);
    */
    
    /* K
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */
    
    /* L
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* M
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */
    
    /* N
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */
    
    /* O
    dlinesegment(uv, vec2(.65,.35), vec2(.65,.65), db);
    da = min(da, db);
    */
    
    /* P
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */
    
    /* Q
    dlinesegment(uv, vec2(-.0,-.65), vec2(.65,-.65), db);
    da = min(da, db);
    */
    
    /* R
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.75), db);
    da = min(da, db);
    */ 
    
    /* S
    dlinesegment(uv, vec2(.65,-.1), vec2(.65,-.65), db);
    da = min(da, db);
    */
    
    /* T
    dlinesegment(uv, vec2(.1,-.75), vec2(.1,0.), db);
    da = min(da, db);
    */
    
    /* U
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* V
    dlinesegment(uv, vec2(-.75,.1), vec2(-.75,.75), db);
    da = min(da, db);
    dlinesegment(uv, vec2(.75,.1), vec2(.75,.75), db);
    da = min(da, db);
    */
    
    /* W
    dlinesegment(uv, vec2(-.65,-.1), vec2(-.65,-.65), db);
    da = min(da, db);
    */
    
    /* X
    dlinesegment(uv, vec2(-.65,-.75), vec2(-.65,-.1), db);
    da = min(da, db);
    */
    
    /* Y
    dlinesegment(uv, vec2(.1,-.75), vec2(.1,-.1), db);
    da = min(da, db);
    */
    
    /* Z
    dlinesegment(uv, vec2(-.65,-.65), vec2(-.65,-.35), db);
    da = min(da, db);
    */
    
    
    
    stroke(da, .1, da);
    smoothmin(d,da,.1,d);

   	col = mix(col, c.xxx, sm(d));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}	

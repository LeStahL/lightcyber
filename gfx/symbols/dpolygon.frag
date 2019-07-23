#version 130
// compute distance to regular polygon
const float pi = acos(-1.);
void dpolygon(in vec2 x, in float N, in float R, out float dst)
{
    float d = 2.*pi/N,
        t = mod(acos(x.x/length(x)), d)-.5*d;
    dst = R-length(x)*cos(t)/cos(.5*d);
}

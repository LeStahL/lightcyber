#version 130
void dcircle(in vec2 x, in float R, out float d)
{
    d = abs(length(x)-R);
}

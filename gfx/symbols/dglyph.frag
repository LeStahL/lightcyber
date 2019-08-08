#version 130

const vec3 c = vec3(1.,0.,-1.);

void rfloat(in float off, out float val);
void dbox(in vec2 x, in vec2 b, out float dst);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
void dcircle(in vec2 x, out float d);
void dcirclesegment(in vec2 x, in float r, in float p0, in float p1, out float d);

void dglyph(in vec2 x, in float ordinal, in float size, out float dst)
{
    float dis;
    dbox(x, 2.*size*c.xx, dis);
    if(dis > 0.)
    {
        dst = dis+.5*size;
        return;
    }

    // Find glyph offset in glyph index
    float nglyphs, offset = 0;
    rfloat(1., nglyphs);
        
    for(float i=0.; i<nglyphs; i+=1.)
    {
        float ord;
        rfloat(2.+2.*i, ord);
        ord = floor(ord);
        
        if(ord == ordinal)
        {
            rfloat(2.+2.*i+1., offset);
            offset = floor(offset);
            break;
        }
    }
    
    if(offset == 0.) 
    {
        dst = 1.;
        return;
    }
    
    // Get distance from glyph data
    float d = 1.;
    
    // Lines
    float nlines;
    rfloat(offset, nlines);
    nlines = floor(nlines);
    offset += 1.;
    for(float i=0.; i<nlines; i+=1.)
    {
        float x1;
        rfloat(offset, x1);
        offset += 1.;
        float y1;
        rfloat(offset, y1);
        offset += 1.;
        float x2;
        rfloat(offset, x2);
        offset += 1.;
        float y2;
        rfloat(offset, y2);
        offset += 1.;
        float da;
        dlinesegment(x, size*vec2(x1,y1), size*vec2(x2, y2), da);
        d = min(d,da);
    }
    
    // Circles
    float ncircles;
    rfloat(offset, ncircles);
    ncircles = floor(ncircles);
    offset += 1.;
    for(float i=0.; i<ncircles; i+=1.)
    {
        float xc;
        rfloat(offset, xc);
        offset += 1.;
        float yc;
        rfloat(offset, yc);
        offset += 1.;
        float r;
        rfloat(offset, r);
        offset += 1.;
        float da;
        dcircle( (x-size*vec2(xc, yc))/size/r,da);
        d = min(d, da*size*r);
    }
    
    // Circle segments
    float nsegments;
    rfloat(offset, nsegments);
    nsegments = floor(nsegments);
    offset += 1.;
    for(float i=0.; i<nsegments; i+=1.)
    {
        float xc;
        rfloat(offset, xc);
        offset += 1.;
        float yc;
        rfloat(offset, yc);
        offset += 1.;
        float r;
        rfloat(offset, r);
        offset += 1.;
        float phi0;
        rfloat(offset, phi0);
        offset += 1.;
        float phi1;
        rfloat(offset, phi1);
        offset += 1.;
        float da;
        dcirclesegment(x-size*vec2(xc,yc), size*r, phi0, phi1, da);
        d = min(d, da);
    }
    
    if(nlines+ncircles+nsegments == 0.)
        dst = dis;
    else dst = d;
}

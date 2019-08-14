#version 130
uniform float iTime;
void scale(out float s)
{
    if(iTime >=  0.0  && iTime <  4.705882 )
    {
        s = mod(iTime- 0.0 , 2.3529 )- 1.17645 ;
        s = smoothstep( -0.196075 ,0.,s)*(1.-smoothstep(0., 0.588225 ,s));
    }
    if(iTime >=  4.705882  && iTime <  18.552036 )
    {
        s = mod(iTime- 4.705882 , 2.3079 )- 1.15395 ;
        s = smoothstep( -0.192325 ,0.,s)*(1.-smoothstep(0., 0.576975 ,s));
    }
    if(iTime >=  18.552036  && iTime <  22.996481 )
    {
        s = mod(iTime- 18.552036 , 2.2222 )- 1.1111 ;
        s = smoothstep( -0.18518333333333334 ,0.,s)*(1.-smoothstep(0., 0.55555 ,s));
    }
    if(iTime >=  22.996481  && iTime <  25.139338 )
    {
        s = mod(iTime- 22.996481 , 2.1427 )- 1.07135 ;
        s = smoothstep( -0.17855833333333335 ,0.,s)*(1.-smoothstep(0., 0.535675 ,s));
    }
    if(iTime >=  25.139338  && iTime <  27.208303 )
    {
        s = mod(iTime- 25.139338 , 2.0691 )- 1.03455 ;
        s = smoothstep( -0.17242500000000002 ,0.,s)*(1.-smoothstep(0., 0.517275 ,s));
    }
    if(iTime >=  27.208303  && iTime <  65.208303 )
    {
        s = mod(iTime- 27.208303 , 2.0 )- 1.0 ;
        s = smoothstep( -0.16666666666666666 ,0.,s)*(1.-smoothstep(0., 0.5 ,s));
    }
    if(iTime >=  65.208303  && iTime <  71.109943 )
    {
        s = mod(iTime- 65.208303 , 1.9673 )- 0.98365 ;
        s = smoothstep( -0.16394166666666668 ,0.,s)*(1.-smoothstep(0., 0.491825 ,s));
    }
    if(iTime >=  71.109943  && iTime <  82.722846 )
    {
        s = mod(iTime- 71.109943 , 1.9354 )- 0.9677 ;
        s = smoothstep( -0.16128333333333333 ,0.,s)*(1.-smoothstep(0., 0.48385 ,s));
    }
    if(iTime >=  82.722846  && iTime <  86.722846 )
    {
        s = mod(iTime- 82.722846 , 2.0 )- 1.0 ;
        s = smoothstep( -0.16666666666666666 ,0.,s)*(1.-smoothstep(0., 0.5 ,s));
    }
    if(iTime >=  86.722846  && iTime <  94.998708 )
    {
        s = mod(iTime- 86.722846 , 2.0691 )- 1.03455 ;
        s = smoothstep( -0.17242500000000002 ,0.,s)*(1.-smoothstep(0., 0.517275 ,s));
    }
    if(iTime >=  94.998708  && iTime <  103.134301 )
    {
        s = mod(iTime- 94.998708 , 2.0338 )- 1.0169 ;
        s = smoothstep( -0.16948333333333332 ,0.,s)*(1.-smoothstep(0., 0.50845 ,s));
    }
    if(iTime >=  103.134301  && iTime <  105.134301 )
    {
        s = mod(iTime- 103.134301 , 2.0 )- 1.0 ;
        s = smoothstep( -0.16666666666666666 ,0.,s)*(1.-smoothstep(0., 0.5 ,s));
    }
    if(iTime >=  105.134301  && iTime <  107.069785 )
    {
        s = mod(iTime- 105.134301 , 1.9354 )- 0.9677 ;
        s = smoothstep( -0.16128333333333333 ,0.,s)*(1.-smoothstep(0., 0.48385 ,s));
    }
    if(iTime >=  107.069785  && iTime <  125.819785 )
    {
        s = mod(iTime- 107.069785 , 1.8751 )- 0.93755 ;
        s = smoothstep( -0.15625833333333333 ,0.,s)*(1.-smoothstep(0., 0.468775 ,s));
    }
    if(iTime >=  125.819785  && iTime <  129.819785 )
    {
        s = mod(iTime- 125.819785 , 2.0 )- 1.0 ;
        s = smoothstep( -0.16666666666666666 ,0.,s)*(1.-smoothstep(0., 0.5 ,s));
    }
    if(iTime >=  129.819785  && iTime <  141.432688 )
    {
        s = mod(iTime- 129.819785 , 1.9354 )- 0.9677 ;
        s = smoothstep( -0.16128333333333333 ,0.,s)*(1.-smoothstep(0., 0.48385 ,s));
    }
    if(iTime >=  141.432688  && iTime <  143.368172 )
    {
        s = mod(iTime- 141.432688 , 1.9354 )- 0.9677 ;
        s = smoothstep( -0.16128333333333333 ,0.,s)*(1.-smoothstep(0., 0.48385 ,s));
    }
    if(iTime >=  143.368172  && iTime <  145.243172 )
    {
        s = mod(iTime- 143.368172 , 1.8751 )- 0.93755 ;
        s = smoothstep( -0.15625833333333333 ,0.,s)*(1.-smoothstep(0., 0.468775 ,s));
    }
    if(iTime >=  145.243172  && iTime <  187.061354 )
    {
        s = mod(iTime- 145.243172 , 1.8182 )- 0.9091 ;
        s = smoothstep( -0.15151666666666666 ,0.,s)*(1.-smoothstep(0., 0.45455 ,s));
    }
}

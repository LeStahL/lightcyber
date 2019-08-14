#version 130
uniform float iTime;
void scale(out float s)
{
    if(iTime >=  0.0  && iTime <  4.705882 )
    {
        s = mod(iTime , 0.588225 )- 0.2941125 ;
        s = smoothstep( -0.04901875 ,0.,s)*(1.-smoothstep(0., 0.14705625 ,s));
    }
    if(iTime >=  4.705882  && iTime <  18.552036 )
    {
        s = mod(iTime , 0.576975 )- 0.2884875 ;
        s = smoothstep( -0.04808125 ,0.,s)*(1.-smoothstep(0., 0.14424375 ,s));
    }
    if(iTime >=  18.552036  && iTime <  22.996481 )
    {
        s = mod(iTime , 0.55555 )- 0.277775 ;
        s = smoothstep( -0.046295833333333335 ,0.,s)*(1.-smoothstep(0., 0.1388875 ,s));
    }
    if(iTime >=  22.996481  && iTime <  25.139338 )
    {
        s = mod(iTime , 0.535675 )- 0.2678375 ;
        s = smoothstep( -0.04463958333333334 ,0.,s)*(1.-smoothstep(0., 0.13391875 ,s));
    }
    if(iTime >=  25.139338  && iTime <  27.208303 )
    {
        s = mod(iTime , 0.517275 )- 0.2586375 ;
        s = smoothstep( -0.043106250000000006 ,0.,s)*(1.-smoothstep(0., 0.12931875 ,s));
    }
    if(iTime >=  27.208303  && iTime <  65.208303 )
    {
        s = mod(iTime , 0.5 )- 0.25 ;
        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));
    }
    if(iTime >=  65.208303  && iTime <  71.109943 )
    {
        s = mod(iTime , 0.491825 )- 0.2459125 ;
        s = smoothstep( -0.04098541666666667 ,0.,s)*(1.-smoothstep(0., 0.12295625 ,s));
    }
    if(iTime >=  71.109943  && iTime <  82.722846 )
    {
        s = mod(iTime , 0.48385 )- 0.241925 ;
        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));
    }
    if(iTime >=  82.722846  && iTime <  86.722846 )
    {
        s = mod(iTime , 0.5 )- 0.25 ;
        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));
    }
    if(iTime >=  86.722846  && iTime <  94.998708 )
    {
        s = mod(iTime , 0.517275 )- 0.2586375 ;
        s = smoothstep( -0.043106250000000006 ,0.,s)*(1.-smoothstep(0., 0.12931875 ,s));
    }
    if(iTime >=  94.998708  && iTime <  103.134301 )
    {
        s = mod(iTime , 0.50845 )- 0.254225 ;
        s = smoothstep( -0.04237083333333333 ,0.,s)*(1.-smoothstep(0., 0.1271125 ,s));
    }
    if(iTime >=  103.134301  && iTime <  105.134301 )
    {
        s = mod(iTime , 0.5 )- 0.25 ;
        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));
    }
    if(iTime >=  105.134301  && iTime <  107.069785 )
    {
        s = mod(iTime , 0.48385 )- 0.241925 ;
        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));
    }
    if(iTime >=  107.069785  && iTime <  125.819785 )
    {
        s = mod(iTime , 0.468775 )- 0.2343875 ;
        s = smoothstep( -0.03906458333333333 ,0.,s)*(1.-smoothstep(0., 0.11719375 ,s));
    }
    if(iTime >=  125.819785  && iTime <  129.819785 )
    {
        s = mod(iTime , 0.5 )- 0.25 ;
        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));
    }
    if(iTime >=  129.819785  && iTime <  141.432688 )
    {
        s = mod(iTime , 0.48385 )- 0.241925 ;
        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));
    }
    if(iTime >=  141.432688  && iTime <  143.368172 )
    {
        s = mod(iTime , 0.48385 )- 0.241925 ;
        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));
    }
    if(iTime >=  143.368172  && iTime <  145.243172 )
    {
        s = mod(iTime , 0.468775 )- 0.2343875 ;
        s = smoothstep( -0.03906458333333333 ,0.,s)*(1.-smoothstep(0., 0.11719375 ,s));
    }
    if(iTime >=  145.243172  && iTime <  187.061354 )
    {
        s = mod(iTime , 0.45455 )- 0.227275 ;
        s = smoothstep( -0.037879166666666665 ,0.,s)*(1.-smoothstep(0., 0.1136375 ,s));
    }
}

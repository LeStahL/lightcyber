#ifndef SCENES_H
#define SCENES_H

    if(t < 11.)
    {
        glUseProgram(logo210_program);
        glUniform1f(logo210_iTime_location, t);
        glUniform2f(logo210_iResolution_location, w, h);
        
#ifdef MIDI
        
        glUniform1f(logo210_iFader0_location, fader0);
        glUniform1f(logo210_iFader1_location, fader1);
        glUniform1f(logo210_iFader2_location, fader2);
        glUniform1f(logo210_iFader3_location, fader3);
        glUniform1f(logo210_iFader4_location, fader4);
        glUniform1f(logo210_iFader5_location, fader5);
        glUniform1f(logo210_iFader6_location, fader6);
        glUniform1f(logo210_iFader7_location, fader7);

        if(override_index == 0)
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 23)
    {
        glUseProgram(evoke_program);
        glUniform1f(evoke_iTime_location, t-11);
        glUniform2f(evoke_iResolution_location, w, h);
        
#ifdef MIDI
        glUniform1f(evoke_iFader0_location, fader0);
        glUniform1f(evoke_iFader1_location, fader1);
        glUniform1f(evoke_iFader2_location, fader2);
        glUniform1f(evoke_iFader3_location, fader3);
        glUniform1f(evoke_iFader4_location, fader4);
        glUniform1f(evoke_iFader5_location, fader5);
        glUniform1f(evoke_iFader6_location, fader6);
        glUniform1f(evoke_iFader7_location, fader7);
        
        if(override_index == 1)
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 39)
    {
        glUseProgram(graffiti_program);
        glUniform1f(graffiti_iTime_location, t-23);
        glUniform2f(graffiti_iResolution_location, w, h);
        
#ifdef MIDI
        glUniform1f(graffiti_iFader0_location, fader0);
        glUniform1f(graffiti_iFader1_location, fader1);
        glUniform1f(graffiti_iFader2_location, fader2);
        glUniform1f(graffiti_iFader3_location, fader3);
        glUniform1f(graffiti_iFader4_location, fader4);
        glUniform1f(graffiti_iFader5_location, fader5);
        glUniform1f(graffiti_iFader6_location, fader6);
        glUniform1f(graffiti_iFader7_location, fader7);
        
        if(override_index == 2)
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 55)
    {
        glUseProgram(groundboxes_program);
        glUniform1f(groundboxes_iTime_location, t-39);
        glUniform2f(groundboxes_iResolution_location, w, h);
        
#ifdef MIDI
        glUniform1f(groundboxes_iFader0_location, fader0);
        glUniform1f(groundboxes_iFader1_location, fader1);
        glUniform1f(groundboxes_iFader2_location, fader2);
        glUniform1f(groundboxes_iFader3_location, fader3);
        glUniform1f(groundboxes_iFader4_location, fader4);
        glUniform1f(groundboxes_iFader5_location, fader5);
        glUniform1f(groundboxes_iFader6_location, fader6);
        glUniform1f(groundboxes_iFader7_location, fader7);
        
        if(override_index == 3) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 71)
    {
        glUseProgram(canal_program);
        glUniform1f(canal_iTime_location, t-55.);
        glUniform2f(canal_iResolution_location, w, h);
        
#ifdef MIDI
        glUniform1f(canal_iFader0_location, fader0);
        glUniform1f(canal_iFader1_location, fader1);
        glUniform1f(canal_iFader2_location, fader2);
        glUniform1f(canal_iFader3_location, fader3);
        glUniform1f(canal_iFader4_location, fader4);
        glUniform1f(canal_iFader5_location, fader5);
        glUniform1f(canal_iFader6_location, fader6);
        glUniform1f(canal_iFader7_location, fader7);
        
        if(override_index == 4)
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 86.5)
    {
        glUseProgram(graffiti_program);
        glUniform1f(graffiti_iTime_location, t);
        glUniform2f(graffiti_iResolution_location, w, h);
        
#ifdef MIDI
        glUniform1f(graffiti_iFader0_location, fader0);
        glUniform1f(graffiti_iFader1_location, fader1);
        glUniform1f(graffiti_iFader2_location, fader2);
        glUniform1f(graffiti_iFader3_location, fader3);
        glUniform1f(graffiti_iFader4_location, fader4);
        glUniform1f(graffiti_iFader5_location, fader5);
        glUniform1f(graffiti_iFader6_location, fader6);
        glUniform1f(graffiti_iFader7_location, fader7);
        
        if(override_index == 5)
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 103)
    {
        glUseProgram(voronoidesign_program);
        glUniform1f(voronoidesign_iTime_location, t-86.5);
        glUniform2f(voronoidesign_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(voronoidesign_iFader0_location, fader0);
        glUniform1f(voronoidesign_iFader1_location, fader1);
        glUniform1f(voronoidesign_iFader2_location, fader2);
        glUniform1f(voronoidesign_iFader3_location, fader3);
        glUniform1f(voronoidesign_iFader4_location, fader4);
        glUniform1f(voronoidesign_iFader5_location, fader5);
        glUniform1f(voronoidesign_iFader6_location, fader6);
        glUniform1f(voronoidesign_iFader7_location, fader7);
        
        if(override_index == 6) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 118.)
    {
        glUseProgram(transbubbles_program);
        glUniform1f(transbubbles_iTime_location, t-103);
        glUniform2f(transbubbles_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(transbubbles_iFader0_location, fader0);
        glUniform1f(transbubbles_iFader1_location, fader1);
        glUniform1f(transbubbles_iFader2_location, fader2);
        glUniform1f(transbubbles_iFader3_location, fader3);
        glUniform1f(transbubbles_iFader4_location, fader4);
        glUniform1f(transbubbles_iFader5_location, fader5);
        glUniform1f(transbubbles_iFader6_location, fader6);
        glUniform1f(transbubbles_iFader7_location, fader7);
        
        if(override_index == 7) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 133.)
    {
        glUseProgram(volclouds_program);
        glUniform1f(volclouds_iTime_location, t-118);
        glUniform2f(volclouds_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(volclouds_iFader0_location, fader0);
        glUniform1f(volclouds_iFader1_location, fader1);
        glUniform1f(volclouds_iFader2_location, fader2);
        glUniform1f(volclouds_iFader3_location, fader3);
        glUniform1f(volclouds_iFader4_location, fader4);
        glUniform1f(volclouds_iFader5_location, fader5);
        glUniform1f(volclouds_iFader6_location, fader6);
        glUniform1f(volclouds_iFader7_location, fader7);
        
        if(override_index == 8) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 141.)
    {
        glUseProgram(chart_program);
        glUniform1f(chart_iTime_location, t-133);
        glUniform2f(chart_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(chart_iFader0_location, fader0);
        glUniform1f(chart_iFader1_location, fader1);
        glUniform1f(chart_iFader2_location, fader2);
        glUniform1f(chart_iFader3_location, fader3);
        glUniform1f(chart_iFader4_location, fader4);
        glUniform1f(chart_iFader5_location, fader5);
        glUniform1f(chart_iFader6_location, fader6);
        glUniform1f(chart_iFader7_location, fader7);
        
        if(override_index == 9) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.)
    {
        glUseProgram(greet_program);
        glUniform1f(greet_iTime_location, t-141);
        glUniform2f(greet_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(greet_iFader0_location, fader0);
        glUniform1f(greet_iFader1_location, fader1);
        glUniform1f(greet_iFader2_location, fader2);
        glUniform1f(greet_iFader3_location, fader3);
        glUniform1f(greet_iFader4_location, fader4);
        glUniform1f(greet_iFader5_location, fader5);
        glUniform1f(greet_iFader6_location, fader6);
        glUniform1f(greet_iFader7_location, fader7);
        
        if(override_index == 10) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    // FIXME hier geballer, taktlänge 1.8182
    else if(t < 156.+1.*1.8182)
    {
        glUseProgram(volclouds_program);
        glUniform1f(volclouds_iTime_location, t);
        glUniform2f(volclouds_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(volclouds_iFader0_location, fader0);
        glUniform1f(volclouds_iFader1_location, fader1);
        glUniform1f(volclouds_iFader2_location, fader2);
        glUniform1f(volclouds_iFader3_location, fader3);
        glUniform1f(volclouds_iFader4_location, fader4);
        glUniform1f(volclouds_iFader5_location, fader5);
        glUniform1f(volclouds_iFader6_location, fader6);
        glUniform1f(volclouds_iFader7_location, fader7);
        
        if(override_index == 11) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+2.*1.8182)
    {
        glUseProgram(graffiti_program);
        glUniform1f(graffiti_iTime_location, t);
        glUniform2f(graffiti_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(graffiti_iFader0_location, fader0);
        glUniform1f(graffiti_iFader1_location, fader1);
        glUniform1f(graffiti_iFader2_location, fader2);
        glUniform1f(graffiti_iFader3_location, fader3);
        glUniform1f(graffiti_iFader4_location, fader4);
        glUniform1f(graffiti_iFader5_location, fader5);
        glUniform1f(graffiti_iFader6_location, fader6);
        glUniform1f(graffiti_iFader7_location, fader7);
        
        if(override_index == 12) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+3.*1.8182)
    {
        glUseProgram(transbubbles_program);
        glUniform1f(transbubbles_iTime_location, t);
        glUniform2f(transbubbles_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(transbubbles_iFader0_location, fader0);
        glUniform1f(transbubbles_iFader1_location, fader1);
        glUniform1f(transbubbles_iFader2_location, fader2);
        glUniform1f(transbubbles_iFader3_location, fader3);
        glUniform1f(transbubbles_iFader4_location, fader4);
        glUniform1f(transbubbles_iFader5_location, fader5);
        glUniform1f(transbubbles_iFader6_location, fader6);
        glUniform1f(transbubbles_iFader7_location, fader7);
        
        if(override_index == 13) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+3.5*1.8182)
    {
        glUseProgram(voronoidesign_program);
        glUniform1f(voronoidesign_iTime_location, t);
        glUniform2f(voronoidesign_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(voronoidesign_iFader0_location, fader0);
        glUniform1f(voronoidesign_iFader1_location, fader1);
        glUniform1f(voronoidesign_iFader2_location, fader2);
        glUniform1f(voronoidesign_iFader3_location, fader3);
        glUniform1f(voronoidesign_iFader4_location, fader4);
        glUniform1f(voronoidesign_iFader5_location, fader5);
        glUniform1f(voronoidesign_iFader6_location, fader6);
        glUniform1f(voronoidesign_iFader7_location, fader7);
        
        if(override_index == 14) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+4.*1.8182)
    {
        glUseProgram(groundboxes_program);
        glUniform1f(groundboxes_iTime_location, t);
        glUniform2f(groundboxes_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(groundboxes_iFader0_location, fader0);
        glUniform1f(groundboxes_iFader1_location, fader1);
        glUniform1f(groundboxes_iFader2_location, fader2);
        glUniform1f(groundboxes_iFader3_location, fader3);
        glUniform1f(groundboxes_iFader4_location, fader4);
        glUniform1f(groundboxes_iFader5_location, fader5);
        glUniform1f(groundboxes_iFader6_location, fader6);
        glUniform1f(groundboxes_iFader7_location, fader7);
        
        if(override_index == 15) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+5.*1.8182)
    {
        glUseProgram(volclouds_program);
        glUniform1f(volclouds_iTime_location, t);
        glUniform2f(volclouds_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(volclouds_iFader0_location, fader0);
        glUniform1f(volclouds_iFader1_location, fader1);
        glUniform1f(volclouds_iFader2_location, fader2);
        glUniform1f(volclouds_iFader3_location, fader3);
        glUniform1f(volclouds_iFader4_location, fader4);
        glUniform1f(volclouds_iFader5_location, fader5);
        glUniform1f(volclouds_iFader6_location, fader6);
        glUniform1f(volclouds_iFader7_location, fader7);
        
        if(override_index == 16) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+6.*1.8182)
    {
        glUseProgram(transbubbles_program);
        glUniform1f(transbubbles_iTime_location, t);
        glUniform2f(transbubbles_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(transbubbles_iFader0_location, fader0);
        glUniform1f(transbubbles_iFader1_location, fader1);
        glUniform1f(transbubbles_iFader2_location, fader2);
        glUniform1f(transbubbles_iFader3_location, fader3);
        glUniform1f(transbubbles_iFader4_location, fader4);
        glUniform1f(transbubbles_iFader5_location, fader5);
        glUniform1f(transbubbles_iFader6_location, fader6);
        glUniform1f(transbubbles_iFader7_location, fader7);
        
        if(override_index == 17) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+8.*1.8182)
    {
        glUseProgram(canal_program);
        glUniform1f(canal_iTime_location, t);
        glUniform2f(canal_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(canal_iFader0_location, fader0);
        glUniform1f(canal_iFader1_location, fader1);
        glUniform1f(canal_iFader2_location, fader2);
        glUniform1f(canal_iFader3_location, fader3);
        glUniform1f(canal_iFader4_location, fader4);
        glUniform1f(canal_iFader5_location, fader5);
        glUniform1f(canal_iFader6_location, fader6);
        glUniform1f(canal_iFader7_location, fader7);
        
        if(override_index == 18) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+9.*1.8182)
    {
        glUseProgram(groundboxes_program);
        glUniform1f(groundboxes_iTime_location, t);
        glUniform2f(groundboxes_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(groundboxes_iFader0_location, fader0);
        glUniform1f(groundboxes_iFader1_location, fader1);
        glUniform1f(groundboxes_iFader2_location, fader2);
        glUniform1f(groundboxes_iFader3_location, fader3);
        glUniform1f(groundboxes_iFader4_location, fader4);
        glUniform1f(groundboxes_iFader5_location, fader5);
        glUniform1f(groundboxes_iFader6_location, fader6);
        glUniform1f(groundboxes_iFader7_location, fader7);
        
        if(override_index == 19) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+10.*1.8182)
    {
        glUseProgram(voronoidesign_program);
        glUniform1f(voronoidesign_iTime_location, t);
        glUniform2f(voronoidesign_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(voronoidesign_iFader0_location, fader0);
        glUniform1f(voronoidesign_iFader1_location, fader1);
        glUniform1f(voronoidesign_iFader2_location, fader2);
        glUniform1f(voronoidesign_iFader3_location, fader3);
        glUniform1f(voronoidesign_iFader4_location, fader4);
        glUniform1f(voronoidesign_iFader5_location, fader5);
        glUniform1f(voronoidesign_iFader6_location, fader6);
        glUniform1f(voronoidesign_iFader7_location, fader7);
        
        if(override_index == 20) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+11.*1.8182)
    {
        glUseProgram(canal_program);
        glUniform1f(canal_iTime_location, t);
        glUniform2f(canal_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(canal_iFader0_location, fader0);
        glUniform1f(canal_iFader1_location, fader1);
        glUniform1f(canal_iFader2_location, fader2);
        glUniform1f(canal_iFader3_location, fader3);
        glUniform1f(canal_iFader4_location, fader4);
        glUniform1f(canal_iFader5_location, fader5);
        glUniform1f(canal_iFader6_location, fader6);
        glUniform1f(canal_iFader7_location, fader7);
        
        if(override_index == 21) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+11.5*1.8182)
    {
        glUseProgram(graffiti_program);
        glUniform1f(graffiti_iTime_location, t);
        glUniform2f(graffiti_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(graffiti_iFader0_location, fader0);
        glUniform1f(graffiti_iFader1_location, fader1);
        glUniform1f(graffiti_iFader2_location, fader2);
        glUniform1f(graffiti_iFader3_location, fader3);
        glUniform1f(graffiti_iFader4_location, fader4);
        glUniform1f(graffiti_iFader5_location, fader5);
        glUniform1f(graffiti_iFader6_location, fader6);
        glUniform1f(graffiti_iFader7_location, fader7);
        
        if(override_index == 22) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+12.*1.8182)
    {
        glUseProgram(transbubbles_program);
        glUniform1f(transbubbles_iTime_location, t);
        glUniform2f(transbubbles_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(transbubbles_iFader0_location, fader0);
        glUniform1f(transbubbles_iFader1_location, fader1);
        glUniform1f(transbubbles_iFader2_location, fader2);
        glUniform1f(transbubbles_iFader3_location, fader3);
        glUniform1f(transbubbles_iFader4_location, fader4);
        glUniform1f(transbubbles_iFader5_location, fader5);
        glUniform1f(transbubbles_iFader6_location, fader6);
        glUniform1f(transbubbles_iFader7_location, fader7);
        
        if(override_index == 23) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+13.*1.8182)
    {
        glUseProgram(volclouds_program);
        glUniform1f(volclouds_iTime_location, t);
        glUniform2f(volclouds_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(volclouds_iFader0_location, fader0);
        glUniform1f(volclouds_iFader1_location, fader1);
        glUniform1f(volclouds_iFader2_location, fader2);
        glUniform1f(volclouds_iFader3_location, fader3);
        glUniform1f(volclouds_iFader4_location, fader4);
        glUniform1f(volclouds_iFader5_location, fader5);
        glUniform1f(volclouds_iFader6_location, fader6);
        glUniform1f(volclouds_iFader7_location, fader7);
        
        if(override_index == 24) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < 156.+14.*1.8182)
    {
        glUseProgram(groundboxes_program);
        glUniform1f(groundboxes_iTime_location, t);
        glUniform2f(groundboxes_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(groundboxes_iFader0_location, fader0);
        glUniform1f(groundboxes_iFader1_location, fader1);
        glUniform1f(groundboxes_iFader2_location, fader2);
        glUniform1f(groundboxes_iFader3_location, fader3);
        glUniform1f(groundboxes_iFader4_location, fader4);
        glUniform1f(groundboxes_iFader5_location, fader5);
        glUniform1f(groundboxes_iFader6_location, fader6);
        glUniform1f(groundboxes_iFader7_location, fader7);
        
        if(override_index == 25) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else if(t < t_end)
    {
        glUseProgram(voronoidesign_program);
        glUniform1f(voronoidesign_iTime_location, t);
        glUniform2f(voronoidesign_iResolution_location, w, h);
      
#ifdef MIDI
        glUniform1f(voronoidesign_iFader0_location, fader0);
        glUniform1f(voronoidesign_iFader1_location, fader1);
        glUniform1f(voronoidesign_iFader2_location, fader2);
        glUniform1f(voronoidesign_iFader3_location, fader3);
        glUniform1f(voronoidesign_iFader4_location, fader4);
        glUniform1f(voronoidesign_iFader5_location, fader5);
        glUniform1f(voronoidesign_iFader6_location, fader6);
        glUniform1f(voronoidesign_iFader7_location, fader7);
        
        if(override_index == 26) 
        {
            select_button(override_index);
            scene_override = 0;
        }
#endif
    }
    else ExitProcess(0);

#endif

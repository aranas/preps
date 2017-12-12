####### INITIATION #######
scenario 		= "preps_exp";
pcl_file 		= "preps_expMAIN.pcl";

active_buttons 	= 4;
# 1 = space
# 2 = enter
# 3 = green / 1 
# 4 = yellow / 0

#---------------------- display specification -------------------------
default_background_color = 0, 0, 0; # black
default_font = "22-Courier";
default_font_size = 25;
default_text_color = 235, 235, 235; # light grey
default_text_align = align_center;


begin;   

#---------------------- screens ----------------------------------------


picture {} default;

picture { text { caption = " "; } t_Info1; x = 0; y = 0;
          text { caption = " "; } t_Info2; x = 0; y = -40;
          text { caption = "Press [ENTER] to confirm or [Esc] to abort. "; } t_info; x = 0; y = -150;
} p_info; # for getting the ppnr

#----------------------Instruction screens ------------------------------------	
picture {
	text {caption = "Wilkommen bei diesem Experiment! 
Ihre Aufgabe ist es, Sätze zu lesen, die Wort für Wort auf dem Bildschirm präsentiert werden. 
Gelegentlich wird es zu den Sätzen inhaltliche Fragen geben. Lesen Sie deshalb aufmerksam und vorsichtig.

Bitte drücken Sie einen der Knöpfe um fortzufahren"; font_size = 25; } t_instruction1;
	x = 0; y = 0;
	}p_instruction1;
	
picture {
	text {caption =
"Das Experiment ist in vier Blöcke unterteilt. 
Zwischen den Blöcken haben Sie die Möglichkeit, eine Pause zu machen. 
Nutzten Sie die Pause, um sich zu entspannen und bestimmen Sie selbst, wann der nächste Block beginnen soll. 

Bitte drücken Sie einen der Knöpfe um fortzufahren."; font_size = 25; } t_instruction2;
	x = 0; y = 0;
	} p_instruction2;

picture {
	text {caption =
"Die Sätze werden Ihnen Wort für Wort auf dem Bildschirm präsentiert. Vor jedem Satz erscheint ein Kreuz,
fixieren Sie dieses bitte. 
Wenn möglich, sollten Sie nicht blinzeln während die Sätze präsentiert werden.
Nach jedem Satz erscheint ein grünes Kreuz. Wenn das grüne Kreuz erscheint, können Sie blinzeln. 
Wechselt das Kreuz zu weiss, dann beginnt in Kürze der nächste Satz. Blinzeln Sie nun bitte nicht mehr!


Bitte drücken Sie einen der Knöpfe um fortzufahren"; font_size = 25; } t_instruction3;
	x = 0; y = 0;
	}p_instruction3;

picture {
	text {caption =
"Nach einem Satz wird ab und an eine inhaltliche Frage gestellt. Mit der Frage werden Ihnen zwei mögliche Antworten präsentiert.
Drücken Sie auf den linken Knopf, um die linke Antwort auszuwählen und auf den rechten Knopf, um die rechte Antwort auszuwählen. Bitte benutzen Sie
dafür den Zeigefinger der jeweiligen Hand.

Wir beginnen nun mit einigen Übungsversuchen. 

Bereit?
Bitte drücken Sie einen der Knöpfe, um mit den Übungsversuchen zu beginnen."; font_size = 25; } t_instruction4;
	x = 0; y = 0;
	}p_instruction4;

picture {
	text {caption =
"Ende der Übungsversuche.
Wir beginnen nun mit dem Experiment. Sollte es noch Unklarheiten geben, wenden Sie sich bitte jetzt an den Versuchsleiter.

Drücken Sie einen der Knöpfe um mit dem Experiment zu beginnen."; font_size = 25; } t_instruction5;
	x = 0; y = 0;
	}p_instruction5;

picture { 
		text {caption = "Ende dieses Blocks. Bitte nehmen Sie sich Zeit für eine Pause!

	Drücken Sie einen der Knöpfe, um mit dem nächsten Block zu beginnen."; font_size = 25; } t_endofblock;
		x =  0; y = 0;
	} p_endofblock;	


picture { 
		text {caption = "Ende des Experiments. Der Versuchsleiter wird nun zu Ihnen kommen.
		
	Vielen Dank für Ihre Teilnahme!!"; font_size = 25; } t_endofexperiment;
		x =  0; y = 0;
	}p_endofexperiment;

picture { 
		text {caption = "Bitte warten Sie..."; font_size = 25; } t_waitmsg;
		x =  0; y = 0;
	}p_waitmsg;




#---------------------- Trial screens ------------------------------------
#White fixation cross
picture {box { height = 1; width = 30; color = 255,255,255; } b_Horz;
			x = 0; y = 0;
			box { height = 30; width = 1; color = 255,255,255; } b_Vert;
			x = 0; y = 0;
} p_fixation;
#Green fixation cross

picture {box { height = 1; width = 30; color = 19,214,113; } b_HorzGreen;
			x = 0; y = 0;
			box { height = 30; width = 1; color = 19,214,113; } b_VertGreen;
			x = 0; y = 0;
} p_blink;

picture {text { caption = "Question"; font_size = 30; } t_question; x = 0; y = 0;
					text { caption = "Answer1"; font_size = 30;} t_answer1; x = -200; y = -80; 
					text { caption = "Answer2"; font_size = 30;} t_answer2; x = 200; y = -80;
} p_question;
  

picture {text { caption = "Word"; font_size = 30; } t_word1; x = 0; y = 0;
} p_word1;
picture {text { caption = "Word"; font_size = 30; } t_word2; x = 0; y = 0;
} p_word2;
picture {text { caption = "Word"; font_size = 30; } t_word3; x = 0; y = 0;
} p_word3;
picture {text { caption = "Word"; font_size = 30; } t_word4; x = 0; y = 0;
} p_word4;
picture {text { caption = "Word"; font_size = 30; } t_word5; x = 0; y = 0;
} p_word5;
picture {text { caption = "Word"; font_size = 30; } t_word6; x = 0; y = 0;
} p_word6;
picture {text { caption = "Word"; font_size = 30; } t_word7; x = 0; y = 0;
} p_word7;
picture {text { caption = "Word"; font_size = 30; } t_word8; x = 0; y = 0;
} p_word8;
picture {text { caption = "Word"; font_size = 30; } t_word9; x = 0; y = 0;
} p_word9;

####### INITIATION #######
scenario 		= "prepositional_exp";
pcl_file 		= "prepositional_expMAIN.pcl";

active_buttons 	= 3;
button_codes 	= 1,2,3;  # give each button as set up in the response panel an integer number . button 3 should correspond to SPACE on keyboard

no_logfile 		= false;

response_logging 	= log_active ;  #only record responses when reponses is expected (i.e. during question trials)
response_matching = simple_matching; 

#---------------------- display specification -------------------------

default_background_color = 0, 0, 0; # black 
default_font = "arial"; #"Helvetica";
default_font_size = 20;
default_text_color = 235, 235, 235; # light grey
default_text_align = align_center;


#---------------------- duration specifications ------------------------
$gap_duration 		= 145; #IWI	#this is only starting value, will be overwritten by expVARS.pcl
$isi_duration 		= 750; #ITI # jitter between 500 - 1000, overwritten in expSUBS.pcl
$blink_duration	= 1000;

begin;   

#---------------------- screens ----------------------------------------


picture {} default;

#----------------------Instructions ------------------------------------
trial {
	trial_duration = forever;
	trial_type = first_response;	
	picture {
	text {caption = "Wilkommen bei diesem Experiment! 
Ihre Aufgabe ist es, Sätze zu lesen, die Wort für Wort auf dem Bildschirm präsentiert werden. 
Gelegentlich wird es zu den Sätzen inhaltliche Fragen geben. Lesen Sie deshalb aufmerksam und vorsichtig.

Bitte drücken Sie einen der Knöpfe um fortzufahren"; font_size = 20; };
	x = 0; y = 0;
	};
	time = 0;
	code = "Instruction1";
	} Instruction1;
	
trial {
	trial_duration = forever;
	trial_type = first_response;
	picture {
	text {caption =
"Das Experiment ist in vier Blöcke unterteilt. 
Zwischen den Blöcken haben Sie die Möglichkeit, eine Pause zu machen. 
Nutzten Sie die Pause, um sich zu entspannen und bestimmen Sie selbst, wann der nächste Block beginnen soll. 

Bitte drücken Sie einen der Knöpfe um fortzufahren."; font_size = 20; };
	x = 0; y = 0;
	};
	time = 0;
	code = "Instruction2";
} Instruction2;

trial {
	trial_duration = forever;
	trial_type = first_response;
	picture {
	text {caption =
"Die Sätze werden Ihnen Wort für Wort auf dem Bildschirm präsentiert. Vor jedem Satz erscheint ein Kreuz,
fixieren Sie dieses bitte. 
Wenn möglich, sollten Sie nicht blinzeln während die Sätze präsentiert werden.
Nach jedem Satz erscheint ein grünes Kreuz. Wenn das grüne Kreuz erscheint, können Sie blinzeln. 
Wechselt das Kreuz zu weiss, dann beginnt in Kürze der nächste Satz. Blinzeln Sie nun bitte nicht mehr!


Bitte drücken Sie einen der Knöpfe um fortzufahren"; font_size = 20; };
	x = 0; y = 0;
	};
	time = 0;
	code = "Instruction3";
} Instruction3;

trial {
	trial_duration = forever;
	trial_type = first_response;
	picture {
	text {caption =
"Nach einem Satz wird ab und an eine inhaltliche Frage gestellt. Mit der Frage werden Ihnen zwei mögliche Antworten präsentiert.
Drücken Sie auf den linken Knopf, um die linke Antwort auszuwählen und auf den rechten Knopf, um die rechte Antwort auszuwählen. Bitte benutzen Sie
dafür den Zeigefinger der jeweiligen Hand.

Wir beginnen nun mit einigen Übungsversuchen. 

Bereit?
Bitte drücken Sie einen der Knöpfe, um mit den Übungsversuchen zu beginnen."; font_size = 20; };
	x = 0; y = 0;
	};
	time = 0;
	code = "Instruction4";
} Instruction4;

trial {
	trial_duration = forever;
	trial_type = first_response;
	picture {
	text {caption =
"Ende der Übungsversuche.
Wir beginnen nun mit dem Experiment. Sollte es noch Unklarheiten geben, wenden Sie sich bitte jetzt an den Versuchsleiter.

Drücken Sie einen der Knöpfe um mit dem Experiment zu beginnen."; font_size = 20; };
	x = 0; y = 0;
	};
	time = 0;
	code = "Instruction5";
} Instruction5;

trial {
   trial_duration = forever;
	trial_type = first_response;  
   picture { 
		text {caption = "Ende dieses Blocks. Bitte nehmen Sie sich Zeit für eine Pause!

	Drücken Sie einen der Knöpfe, um mit dem nächsten Block zu beginnen."; font_size = 20; } ;
		x =  0; y = 0;
	};	
   time = 0;
} Trial_End_block;

trial {
   trial_duration = forever;
	trial_type = first_response;  
   picture { 
		text {caption = "Ende des Experiments. Der Versuchsleiter wird nun zu Ihnen kommen.
		
	Vielen Dank für Ihre Teilnahme!!"; font_size = 20; } ;
		x =  0; y = 0;
	};	
   time = 0;
} End_experiment;

trial {
   trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 3; 
   picture { 
		text {caption = "Bitte warten Sie..."; font_size = 20; } ;
		x =  0; y = 0;
	};	
   time = 0;
	code = "Waiting_screen";
} Wait_screen;




#----------------------Display specifications for trials ------------------------------------
trial { #showing fixation cross
	all_responses = false;
	stimulus_event {
		picture {
			text { caption = "+"; font_size = 26; font = "Helvetica"; };  
			x = 0; y = 0;
		};
		duration = $isi_duration;
	code = "fixation";		
	} Event_iti;	
} Trial_iti_fix;

trial { #green cross
	all_responses = false;
	stimulus_event {
		picture {
			text { caption = "+"; font_size = 26; font = "Helvetica"; font_color = 19,214,113; };  
			x = 0; y = 0;
		};
		duration = $blink_duration;
	code = "fixation_green";		
	} Event_blink;	
} Trial_iti_blink;


trial { #showing sentence word by word
	all_responses = false;
   start_delay = 0;
   stimulus_event {
		picture {
			text {caption = "xx"; font_size = 26; } Text_words; 
			x = 0; y = 0;
		} Text_sentence; 			 
	} Event_sentence;
   stimulus_event {
		picture default; 			
		duration = $gap_duration; 
   } Sentence_isi;  
} Trial_show_sentence;

trial { #showing question
	start_delay = 0;
	trial_duration = forever;
	trial_type = first_response;
   stimulus_event {
		picture {
			text {caption = "xx"; font_size = 26; } Question_words;
				x = 0; y = 0;
			text { caption = "First word"; font_size = 20;} Answer_1; x = -200; y = -80;
			text { caption = "Second word"; font_size = 20;} Answer_2; x = 200; y = -80;
			} Text_Question; 
		stimulus_time_in = 0;
		stimulus_time_out = never;
   } Question;  
} Trial_show_Question;  


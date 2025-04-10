globals
[
  ; The length of the gene
  buffer ; int
  start_site ; int
  end_site ; int

  ; The state of the promoter (0, ..., 6)
  promoter_state ; int

  ; Size of a polymerase(in nt)
  L ; int

  ; RNA count
  rna_count

  ; Misc. variables
  time ; float
  min_dist ; float
  i ; int
  file_name ; str
  update ; float
  n ; int
  closest ; polymerase
  effect ; float
  cur_xcor ; coord
  topo_xcor ; coord
  x ; float
  y ; float
]

breed [polymerases polymerase]
breed [PICs PIC]
polymerases-own [
  ; Phase of elongation (-4 -3 -2 -1 0 1 2 3 4)
  pol_phase ; int

  ; Amount of supercoiling behind the polymerase
  pol_supercoil_us ; float

  ; Amount of supercoiling ahead of the polymerase
  pol_supercoil_ds ; float
]



to setup
  clear-all
  reset-ticks
  file-close-all


  ; Simulation settings
  set buffer 20
  set start_site buffer
  set end_site (max-pxcor - buffer)
  set promoter_state 0
  set L 15
  set-default-shape polymerases "circle"
  set-default-shape PICs "circle"
  set time 0.0
  set i 0
  set update 500
  set n 1000

  ; Draw gene
  ask patches with [buffer <= pxcor and 19 <= pycor and pycor <= 21 and pxcor <= end_site] [
    set pcolor white
  ]

  ; Draw promoter region (TSS)
  ask patches with [19 <= pycor and pycor <= 21 and pxcor = start_site] [
    set pcolor green
  ]

  ; Draw termintation site
  ask patches with [19 <= pycor and pycor <= 21 and pxcor = end_site] [
    set pcolor red
  ]

  open_file
end

to go
  tick
  set time ticks * dt

  ; Compute minimum polymerase distance(to promoter)
  if (count polymerases) != 0 [
    set min_dist min [xcor] of polymerases
  ]

  ; Transcription Initiation
  if (count polymerases) = 0 or (L + start_site) <= min_dist [ ; Ensure promoter region is clear
    initiate
  ]

  ; Transcription Elongation
  ask polymerases [
    elongate
  ]

  ; Death Process for mRNAs
  death_process

  ; Cancelling out strain
  cancel_strain

  ; Simulate effect of topoisomerase
  topoisomerase


  ; Sampling from Stationary Distribution
  if 5000 <= time and (time) mod update < dt * 0.9 [
    write_out

    ; Print percent progress
    show precision ((time / (5000 + update * n) ) * 100) 1
  ]
  if 5000 + update * n <= time [
    file-close-all
    stop
  ]
end


;; Our Procedures ----
; Default attributes for a polymerase
to def_polymerase
    set size L
    set color blue
    set heading 90
    set pol_phase 0
end

; Transcription Initation
to initiate
  (ifelse
    ; PIC Assembly
    promoter_state = 0 [ ; Unbound promoter
      ; forward | TF2D binds
      if random-float 1.0 < (tf2d_conc * tf2d_kon) * dt [
        ask patch buffer 20 [
          sprout-PICs 1 [
            set size L - 6
            set color yellow
            set heading 90
          ]
        ]
        set promoter_state 0 + 1
      ]
    ]
    promoter_state = 1 [ ; TF2D bound promoter
      ; forward | TF2A binds
      if random-float 1.0 < (tf2a_conc * tf2a_kon) * dt [
        ask PICs [
          set size L - 5
          set color green
        ]
        set promoter_state 1 + 1
      ]

      ; reverse | TF2D dissociates
      if promoter_state = 1 and random-float 1.0 < (tf2d_koff) * dt [
        ask PICs [
          set size L - 6
          set color yellow
        ]
        set promoter_state 1 - 1
      ]
    ]
    promoter_state = 2 [ ; TF2D, TF2A bound promoter
      ; forward | TF2B binds
      if random-float 1.0 < (tf2b_conc * tf2b_kon) * dt [
        ask PICs [
          set size L - 4
          set color lime
        ]
        set promoter_state 2 + 1
      ]

      ; reverse | TF2A dissociates
      if promoter_state = 2 and random-float 1.0 < (tf2a_koff) * dt [
        ask PICs [
          set size L - 5
          set color green
        ]
        set promoter_state 2 - 1]
    ]
    promoter_state = 3 [ ; TF2D, TF2A, TF2B bound promoter
      ; forward | RNAP-TF2F binds
      if random-float 1.0 < (rnap_tf2f_conc * rnap_tf2f_kon) * dt [
        ask PICs [
          set size L - 3
          set color turquoise
        ]

        set promoter_state 3 + 1]

      ; reverse | TF2B dissociates
      if promoter_state = 3 and random-float 1.0 < (tf2b_koff) * dt [
        ask PICs [
          set size L - 4
          set color lime
        ]
        set promoter_state 3 - 1]
    ]
    promoter_state = 4 [ ; TF2D, TF2A, TF2B, RNAP-TF2F bound promoter
      ; forward | TF2E binds
      if random-float 1.0 < (tf2e_conc * tf2e_kon) * dt [
        ask PICs [
          set size L - 2
          set color cyan
        ]
        set promoter_state 4 + 1
      ]

      ; reverse | RNAP-TF2F dissociates
      if promoter_state = 4 and random-float 1.0 < (rnap_tf2f_koff) * dt [
        ask PICs [
          set size L - 3
          set color turquoise
        ]
        set promoter_state 4 - 1
      ]
    ]
    promoter_state = 5 [ ; TF2D, TF2A, TF2B, RNAP-TF2F, TF2E bound promoter
      ; forward | TF2H binds
      if random-float 1.0 < (tf2h_conc * tf2h_kon) * dt [
        ask PICs [
          set size L - 1
          set color sky
        ]
        set promoter_state 5 + 1
      ]

      ; reverse | TF2E dissociates
      if promoter_state = 5 and random-float 1.0 < (tf2e_koff) * dt [
        ask PICs [
          set size L - 2
          set color turquoise
        ]
        set promoter_state 5 - 1
      ]
    ]
    promoter_state = 6 [ ; TF2D, TF2A, TF2B, RNAP-TF2F, TF2E, TF2H bound promoter
      ; forward | TF2H phosphorylates (ligand)
      if random-float 1.0 < (tf2h_phospho) * dt [
        ask PICs [
          set size L
          set color blue
        ]
        set promoter_state 6 + 1
      ]

      ; reverse | TF2H dissociates
      if promoter_state = 6 and random-float 1.0 < (tf2h_koff) * dt [
        ask PICs [
          set size L - 1
          set color cyan
        ]
        set promoter_state 6 - 1
      ]
    ]
    ; Promoter Escape
    promoter_state = 7 [ ; TF2H-phosphorylated CTD
      ; forward | exchange of GTFs for ETFs (successful promoter escape)
      if random-float 1.0 < (exchange) * dt [
        ; Everything unbinds from promoter
        set promoter_state 0
        ask PICs [
          die
        ]

        ; Spawn new `polymerase` agent
        ask patch start_site 20 [
          sprout-polymerases 1 [
            def_polymerase
          ]
        ]
      ]
    ]


  )
end


; Transcription Elongation
to elongate
  if pol_phase = -4 [ set color brown set label "backtrack" ]
  if pol_phase = -3 [ set color red set label "frayed" ]
  if pol_phase = -2 [ set color orange set label "p-preT" ]
  if pol_phase = -1 [ set color pink set label "ePEC" ]
  if pol_phase =  0 [ set color blue set label "preT" ]
  if pol_phase =  1 [ set color cyan set label "postT" ]
  if pol_phase =  2 [ set color green set label "NTP-bound" ]
  if pol_phase =  3 [ set color lime set label "NTP-aligned" ]
  if pol_phase =  4 [ set color violet set label "catalysis" ]

  (ifelse
    ; On-Pathway Class
    pol_phase = 0 [ ; Pre-translocated State
      set color blue
      ; Attempt to move to post-translocated state if no steric hindrance
      if not (any? (other polymerases in-cone L 10)) [
        ; Calculate effect of torsional strain
        set effect 1 / (1 + exp((pol_supercoil_ds - pol_supercoil_us) - 100))

        if random-float 1 < (translocate_kfor * effect) * dt [
          ; Move forward along DNA strand
          set pol_phase 0 + 1
          forward 1

          ; Transcription Termination
          if xcor >= end_site [
            set rna_count rna_count + 1
            die
          ]

          ; Generate torsional strain
          set pol_supercoil_ds pol_supercoil_ds + 1
          set pol_supercoil_us pol_supercoil_us - 1
        ]
      ]

      ; Move to ePEC state(Paused Polymerase)
      if pol_phase = 0 and random-float 1 < (pausing_kfor) * dt [set pol_phase -1]
    ]
    pol_phase = 1 [ ; Post-Translocated State
      ; Bind correct nucleotide
      if random-float 1 < (binding_kfor * ntp_conc) * dt [set pol_phase 1 + 1]

      ; Move back to pre-translocated state
      if pol_phase = 1 and random-float 1 < (translocate_krev) * dt [
        set pol_phase 1 - 1
        back 1

        ; Undo torsional strain
        set pol_supercoil_ds pol_supercoil_ds - 1
        set pol_supercoil_us pol_supercoil_us + 1
      ]
    ]
    pol_phase = 2 [ ; NTP-bound State
      ; Align nucleotide
      if random-float 1 < (alignment_kfor) * dt [set pol_phase 2 + 1]

      ; Nucleotide dissociates
      if pol_phase = 2 and random-float 1 < (binding_krev) * dt [set pol_phase 2 - 1]
    ]
    pol_phase = 3 [ ; NTP-aligned State
      ; Catalysis occurs
      if random-float 1 < (catalysis_kfor) * dt [set pol_phase 3 + 1]

      ; NTP unaligns
      if pol_phase = 3 and random-float 1 < (alignment_krev) * dt [set pol_phase 3 - 1]
    ]
    pol_phase = 4 [ ; NTP-aligned State
      ; PPi release -> pre-translocated state
      if random-float 1 < (ppi_release) * dt [set pol_phase 0]

      ; Catalysis reverses
      if pol_phase = 4 and random-float 1 < (catalysis_krev) * dt [set pol_phase 4 - 1]
    ]
    ; Offline Class ----
    pol_phase = -1 [ ; ePEC State
      set color red
      ; Backtrack to paused pre-translocated state
      if random-float 1 < (p_halftranslocate_kfor) * dt [set pol_phase -1 - 1]

      ; Unpause to moving pre-translocated state
      if pol_phase = -1 and random-float 1 < (pausing_krev) * dt [set pol_phase 0]
    ]
    pol_phase = -2 [ ; Paused Pre-Translocated State
      ; Move to frayed paused state
      if random-float 1 < (p_fray_kfor) * dt [set pol_phase -2 - 1]

      ; Move back to ePEC state
      if pol_phase = -2 and random-float 1 < (p_halftranslocate_krev) * dt [set pol_phase -2 + 1]
    ]
    pol_phase = -3 [ ; Paused Frayed State
      if random-float 1 < (p_backtrack_kfor) * dt [set pol_phase -3 - 1]

      ; Move to pre-translocated state
      if pol_phase = -3 and random-float 1 < (p_fray_krev) * dt [set pol_phase -3 + 1]
    ]
    pol_phase = -4 [ ; Backtracked State
      ; Move to frayed state
      if random-float 1 < (p_backtrack_krev) * dt [set pol_phase -3 + 1]
    ]
  )
end

; Death of mRNA
to death_process
  if random-float 1 < (death_rate * rna_count) * dt [
    set rna_count rna_count - 1
  ]
end

; Cancelling out torsional strain
to cancel_strain
  ask polymerases [
    ; Find closest polymerases ahead ----
    set cur_xcor xcor
    set closest min-one-of other polymerases with [cur_xcor < xcor] [distance myself]

    if closest != nobody [
      set x pol_supercoil_ds
      ask closest [
        set y pol_supercoil_us
      ]

      ; Cancel out strain
      set pol_supercoil_ds (x + y) / 2.0
      ask closest [
        set pol_supercoil_us (x + y) / 2.0
      ]
    ]



    ; Find closest polymerases behind ----
    set closest min-one-of other polymerases with [xcor < cur_xcor] [distance myself]

    if closest != nobody [
      set x pol_supercoil_us
      ask closest [
        set y pol_supercoil_ds
      ]

      ; Cancel out strain
      set pol_supercoil_us (y + x) / 2.0
      ask closest [
        set pol_supercoil_ds (y + x) / 2.0
      ]
    ]
  ]
end

to topoisomerase
  if random-float 1 < (topo_conc) * dt [
    ; Select random position to bind to
    set topo_xcor random (end_site - start_site) + start_site

    ; Remove torsional strain ahead
    set closest min-one-of polymerases with [topo_xcor < xcor and (xcor - topo_xcor) < 50] [xcor - topo_xcor]
    if closest != nobody [
      ask closest [set pol_supercoil_us 0.0]
    ]

    ; Remove torsional strain behind
    set closest min-one-of polymerases with [xcor < topo_xcor and (topo_xcor - xcor) < 50] [topo_xcor - xcor]
    if closest != nobody [
      ask closest [set pol_supercoil_ds 0.0]
    ]
  ]
end

; Open up next simulation's file
to open_file
  while [file-exists? (word "sim-" i ".txt")] [
    set i i + 1
  ]

  set file_name (word "sim-" i ".txt")
  print file_name

  file-open file_name
end

; Write out mRNA count
to write_out
  file-write rna_count
end
@#$#@#$#@
GRAPHICS-WINDOW
28
518
1386
637
-1
-1
2.7
1
10
1
1
1
0
1
1
1
0
499
0
40
0
0
1
ticks
30.0

BUTTON
13
469
79
506
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
177
472
319
505
dt
dt
0.001
0.05
0.05
0.001
1
NIL
HORIZONTAL

TEXTBOX
279
10
478
28
Transcription Initiation Parameters
12
0.0
0

TEXTBOX
914
10
1124
40
Transcription Elongation Parameters
12
0.0
1

BUTTON
94
470
160
507
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
729
39
848
57
Translocation Step (1)
11
0.0
1

SLIDER
702
57
874
90
translocate_kfor
translocate_kfor
0
10
7.6
0.1
1
NIL
HORIZONTAL

SLIDER
702
98
874
131
translocate_krev
translocate_krev
0
10
5.6
0.1
1
NIL
HORIZONTAL

TEXTBOX
926
38
1013
56
NTP Binding (2)
11
0.0
1

SLIDER
883
57
1055
90
binding_kfor
binding_kfor
0
10
5.3
0.1
1
NIL
HORIZONTAL

SLIDER
883
98
1055
131
binding_krev
binding_krev
0
10
5.2
0.1
1
NIL
HORIZONTAL

TEXTBOX
1104
39
1202
57
NTP Alignment (3)
11
0.0
1

SLIDER
1063
57
1235
90
alignment_kfor
alignment_kfor
0
10
3.2
0.1
1
NIL
HORIZONTAL

SLIDER
1063
98
1235
131
alignment_krev
alignment_krev
0
10
3.5
0.1
1
NIL
HORIZONTAL

TEXTBOX
1293
40
1361
58
Catalysis (4)
11
0.0
1

SLIDER
1242
57
1414
90
catalysis_kfor
catalysis_kfor
0
10
4.1
0.1
1
NIL
HORIZONTAL

SLIDER
1243
98
1415
131
catalysis_krev
catalysis_krev
0
10
4.5
0.1
1
NIL
HORIZONTAL

SLIDER
975
141
1147
174
ppi_release
ppi_release
0
10
10.0
0.1
1
NIL
HORIZONTAL

SLIDER
17
59
165
92
tf2a_conc
tf2a_conc
0
10
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
162
34
360
52
Transcription Factor Concentrations
11
0.0
1

SLIDER
177
59
325
92
tf2b_conc
tf2b_conc
0
10
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
328
156
480
174
PIC Assembly Intrinsic Rates
11
0.0
1

SLIDER
171
102
329
135
rnap_tf2f_conc
rnap_tf2f_conc
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
335
103
492
136
tf2h_conc
tf2h_conc
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
15
101
164
134
tf2e_conc
tf2e_conc
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
334
58
490
91
tf2d_conc
tf2d_conc
0
10
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
31
203
122
221
TF2D Binding (1)
11
0.0
1

TEXTBOX
227
203
319
221
TF2A Binding (2)
11
0.0
1

TEXTBOX
427
204
515
222
TF2B Binding (3)
11
0.0
1

TEXTBOX
27
339
151
357
RNAP-TF2F Binding (4)
11
0.0
1

TEXTBOX
227
337
316
355
TF2E Binding (5)
11
0.0
1

TEXTBOX
407
338
497
356
TF2H Binding (6)
11
0.0
1

SLIDER
0
269
172
302
tf2d_koff
tf2d_koff
0
10
1.5
0.1
1
NIL
HORIZONTAL

SLIDER
0
227
172
260
tf2d_kon
tf2d_kon
0
10
3.0
0.1
1
NIL
HORIZONTAL

SLIDER
183
227
355
260
tf2a_kon
tf2a_kon
0
10
3.0
0.1
1
NIL
HORIZONTAL

SLIDER
184
270
356
303
tf2a_koff
tf2a_koff
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
367
227
539
260
tf2b_kon
tf2b_kon
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
368
271
540
304
tf2b_koff
tf2b_koff
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
0
362
172
395
rnap_tf2f_kon
rnap_tf2f_kon
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
1
405
173
438
rnap_tf2f_koff
rnap_tf2f_koff
0
10
4.0
0.1
1
NIL
HORIZONTAL

SLIDER
184
362
356
395
tf2e_kon
tf2e_kon
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
184
406
356
439
tf2e_koff
tf2e_koff
0
10
3.2
0.1
1
NIL
HORIZONTAL

SLIDER
366
362
538
395
tf2h_kon
tf2h_kon
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
366
407
538
440
tf2h_koff
tf2h_koff
0
10
2.4
0.1
1
NIL
HORIZONTAL

TEXTBOX
601
340
693
358
Promoter Escape
11
0.0
1

SLIDER
559
361
731
394
tf2h_phospho
tf2h_phospho
0
10
0.5
0.1
1
NIL
HORIZONTAL

SLIDER
560
407
732
440
exchange
exchange
0
10
0.5
0.1
1
NIL
HORIZONTAL

TEXTBOX
922
397
1234
415
(Effective) Concentration of Elongation-regulating Factors
11
0.0
1

SLIDER
809
419
981
452
topo_conc
topo_conc
0
10
4.0
0.1
1
NIL
HORIZONTAL

SLIDER
989
419
1161
452
ntp_conc
ntp_conc
0
10
7.6
0.1
1
NIL
HORIZONTAL

TEXTBOX
730
206
824
224
Pausing Step (-1)
11
0.0
1

SLIDER
691
229
863
262
pausing_kfor
pausing_kfor
0
10
8.4
0.1
1
NIL
HORIZONTAL

SLIDER
691
269
863
302
pausing_krev
pausing_krev
0
10
4.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
952
188
1165
206
Paused Elongation Complex Parameters
11
0.0
1

TEXTBOX
906
208
1056
226
Pre-Translocation Step (-2)
11
0.0
1

SLIDER
872
229
1065
262
p_halftranslocate_kfor
p_halftranslocate_kfor
0
10
2.6
0.1
1
NIL
HORIZONTAL

SLIDER
871
270
1065
303
p_halftranslocate_krev
p_halftranslocate_krev
0
10
3.2
0.1
1
NIL
HORIZONTAL

SLIDER
1076
228
1248
261
p_fray_kfor
p_fray_kfor
0
10
1.3
0.1
1
NIL
HORIZONTAL

SLIDER
1077
272
1249
305
p_fray_krev
p_fray_krev
0
10
3.7
0.1
1
NIL
HORIZONTAL

TEXTBOX
1120
209
1212
227
Fraying Step (-3)
11
0.0
1

TEXTBOX
1289
209
1439
227
Backtracking Step (-4)
11
0.0
1

SLIDER
1257
228
1429
261
p_backtrack_kfor
p_backtrack_kfor
0
10
0.6
0.1
1
NIL
HORIZONTAL

SLIDER
1257
272
1429
305
p_backtrack_krev
p_backtrack_krev
0
10
3.8
0.1
1
NIL
HORIZONTAL

SLIDER
1211
474
1383
507
death_rate
death_rate
0
0.1
1.0E-4
0.00001
1
NIL
HORIZONTAL

PLOT
28
651
1011
984
RNA Count
time
mrna_count
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"mrna" 1.0 0 -16777216 true "" "plotxy time rna_count"

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@

   M S B A M L   `   `   `  ���� C  ?UAT_EFDC, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null#clr-namespace:UAT_EFDCUAT_EFDC      UAT_EFDC.Window2   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     @d2http://schemas.microsoft.com/expression/blend/2008     Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     (localclr-namespace:UAT_EFDC       � Title$  Edit Getefdc��$	��600��$	��800��.��� .+��� $
��Left=�$��23��$��
10,217,0,0q� �� TextWrapping$
 Wrap=�$��t��$	��Top=�$	��120���� $
��Left=�$��	0,130,0,0q�$
��Wrap=�$	��Top=�$	��431��$	��526���� $S�M*****************************************************************************���� �� $�** OPTIONS FOR OUTPUT:���� �� $R�L** LAYK        = K>0: DATA AT LAYER NUMBER K TO BE EXPORTED AT TIME=JULTIME ���� �� $;�5**                 0: DEPTH-AVERAGED DATA IS EXPORTED���� �� $?�9**                -1: Get High Frequency output FOR CELLS���� �� $R�L**                -2: Extract data for Time series at a height above bed (m)���� �� $;�5**                -3: LOAD TMP.DAT AND EXPORT TECPLOT���� �� $�**���� �� $H�B** ZOPT          = 1: FOR THE DEPTH UNDER WATER SURFACE IF LAYK=-2���� �� $B�<**                 2: FOR THE HEIGHT ABOVE BOTTOM IF LAYK=-2���� �� $A�;** NDRIFTER         : N1:N2 A SET OF DRIFTER TO GET (X,Y,Z)���� �� $�**���� �� $�**���� �� $:�4** JULTIME          : JULIAN TIME FOR SELECTED LAYER���� �� $8�2**                  > MAXTIME THEN JULTIME=MAXTIME���� �� $3�-**                    0 DATA FOR ALL SNAPSHOT���� �� $A�;** NLOC             : NUMBER OF CELLS TO EXTRACT TIMESERIES���� �� $T�N** ROTA          = 1: (U,V) AT CELL CENTER ROTATED TO TRUE EAST AND TRUE NORTH���� �� $A�;**                 0: (U,V) AT CELL  FACES WITHOUT ROTATION���� �� $4�.** INDEX         = 1: (I,J) OF CELLS ARE GIVEN���� �� $4�.**                 0: (X,Y) OF CELLS ARE GIVEN���� �� $3�-** VPROF         = 1: EXPORT VERTICAL PROFILE���� �� $?�9**               = 0: NO EXPORTATION FOR VERTICAL PROFILE���� �� $�**���� �� $3�-** TECPLOT       = 1: EXPORT DATA FOR TECPLOT���� �� $2�,**               = 0: NO TECPLOT EXPORTATION���� �� $�**���� �� $S�M*****************************************************************************���� $
��Left=�$��	35,46,0,0q�$
��Wrap=�$O��ILAYK    JULTIME      NLOC    ROTA  INDEX  VPROF  TECPLOT   ZOPT  NDRIFTER��$	��Top=�$	��623��$��16(��� -   $	��t10��$
��Left=�$	��431��$��505,130,0,0q�$
 Wrap=�$����15              5            0.1
	15              6            0.1
	15              7            0.1
	15              9            0.1
                 15             8            0.1��$	��Top=�$	��288��$��16(� �� AcceptsReturn .  �� 
AcceptsTab . �� -   $��t1��$
��Left=�$��30��$��	35,71,0,0q�$
 Wrap=�$��1��$	��Top=�$��30��$��16(��� -   $��t2��$
��Left=�$��30��$��
102,71,0,0q�$
 Wrap=�$	��1.0��$	��Top=�$��30��$��16(��� -   $��t3��$
��Left=�$��30��$��
189,71,0,0q�$
 Wrap=�$��5��$	��Top=�$��30��$��16(��� -   $��t4��$
��Left=�$��30��$��
259,71,0,0q�$
 Wrap=�$��1��$	��Top=�$��30��$��16(��� -   $��t5��$
��Left=�$��30��$��
319,71,0,0q�$
 Wrap=�$��1��$	��Top=�$��30��$��16(��� -   $��t6��$
��Left=�$��30��$��
372,71,0,0q�$
 Wrap=�$��0��$	��Top=�$��30��$��16(��� -   $��t7��$
��Left=�$��30��$��
434,71,0,0q�$
 Wrap=�$��0��$	��Top=�$��30��$��16(��� -	   $��t8��$
��Left=�$��30��$��
498,71,0,0q�$
 Wrap=�$��1��$	��Top=�$��30��$��16(��� -
   $��t9��$
��Left=�$��30��$��
573,71,0,0q�$
 Wrap=�$��0��$	��Top=�$��30��$��16(��� $
��Left=�$��531,105,0,0q�$
��Wrap=�$(��"I/X           J/Y          ZINT(m)��$	��Top=�$	��246��$��16(�$��25���� -   $��p��$
��Left=�$��30��$��	102,2,0,0q�$
 Wrap=�$	��Top=�$	��579��$��16(��� $
��Left=�$��35,4,0,0q�$
��Wrap=�$��Path:��$	��Top=�$��20(�
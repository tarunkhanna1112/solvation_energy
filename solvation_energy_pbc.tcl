proc water_scale {} {

	# DETERMING THE UNIQUE WATER MOLECULE AROUND EACH ATOM BASED ON CHARGE WITH THE MINIMUM NUMBER OF WATER AROUND ANY ATOM EQUAL 1

	set f [open "[input1]" "r"]
	set data [read $f]
	close $f

	set g [open "scale" "w"]

	set k 0 

	set n_atoms [input2]

	while { [lindex $data $k] != "@<TRIPOS>ATOM" } { 
		incr k
	}
	incr k
	
	set cmax 0.0
	for {set i 0} {$i < $n_atoms} {incr i} {
		set c($i) [lindex $data [expr { $k + 8 }]]
		if { [expr { abs($c($i)) }] > $cmax } {
			set cmax [expr { abs($c($i)) }]
		}
		incr k 9
	}
	for { set i 0} {$i < $n_atoms} {incr i} {
		set nwat [expr { abs(($c($i) / $cmax)) * [input11] }]
		set nwat [expr { round($nwat) }]
		if { $nwat == 0 } {
			set nwat 1
		}
		puts $g "$nwat"
	}
	close $g
}	

proc nearest_wat {} {

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
	set p1(1) "    "
	set p1(2) "   "
	set p1(3) "  "
	set p1(4) " "
	set p1(5) ""

	set c(4) "   "
	set c(5) "  "
	set c(6) " "
	set c(7) ""

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) ""
	
	set g [open "scale" "r"]
	set data1 [read $g]
	close $g

	set mm [open "Amber_energy" "w"]
	puts $mm "# FRAME WAT PL9 WAT+PL9 HYDRATION_ENERGY"
	set qm [open "Gaussian_energy" "w"]
	puts $qm "# FRAME WAT PL9 WAT+PL9 HYDRATION_ENERGY"

	set start_frame [input3]
	set end_frame [input4]
	set step [input5]
	set n_frames [expr { (($end_frame - $start_frame)/$step) + 1 }]

	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {

		puts ""
		puts "			**** FRAME $fr ****"
		puts ""

		set h1 [open "dummy" "w"]

		set h [open "nwat.pdb" "w"]

		set h2 [open "pl9_only.pdb" "w"]

		# EXECUTING CPPTRAJ

		# FORMING THE INITIAL PDB FILE 

		set inp [open "ai.dat" "w"]
		puts $inp "trajin [input6] $fr $fr"
		puts $inp "trajout frame.rst rst"
		puts $inp "trajout frame.pdb pdb"
		puts $inp "go"
		close $inp

		set k2 0

		set inp [open "ai.dat" "r"]

		exec cpptraj -p [input7] -i ai.dat 

		close $inp

		set rs [open "frame.rst" "r"]
		set rst [read $rs]
		close $rs

		while { $k2 < [llength $rst] } {
			incr k2
		}

		set box_x [lindex $rst [expr { $k2 - 3 }]]
		set box_y [lindex $rst [expr { $k2 - 2 }]]
		set box_z [lindex $rst [expr { $k2 - 1 }]]

		

		set f [open "frame.pdb" "r"]
		set data [read $f]
		close $f

		set k2 0	
		while {[lindex $data $k2] != "TER" || [lindex $data [expr { $k2 - 1 }]] != "H" } {
			if { [lindex $data $k2] == "ATOM" &&  [lindex $data [expr { $k2 + 3 }]] == "MOL" } { 
				set coordx [lindex $data [expr { $k2 + 5 }]] 
				set coordy [lindex $data [expr { $k2 + 6 }]] 
				set coordz [lindex $data [expr { $k2 + 7 }]]

				puts $h1 "{ $coordx $coordy $coordz }"

				set rn1 [lindex $data [expr { $k2 + 1 }]]
				set srn1 [string length $rn1]

				set at1 [lindex $data [expr { $k2 + 2 }]]
				set sat1 [string length $at1]

				set an1 [lindex $data [expr { $k2 + 4 }]]
				set san1 [string length $an1]

				set x1 [lindex $data [expr { $k2 + 5 }]]
				set sx1 [string length $x1]

				set y1 [lindex $data [expr { $k2 + 6 }]] 
				set sy1 [string length $y1]

				set z1 [lindex $data [expr { $k2 + 7 }]] 
				set sz1 [string length $z1]

				puts $h2 "ATOM   $p($srn1)$rn1 $at1$sat($sat1) [input19]  $p($san1)$an1     $c($sx1)$x1 $c($sy1)$y1 $c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k2 + 10 }]]"

			}
			incr k2
		}
		close $h1
		puts $h2 "TER"
		close $h2

		# CALCULATING THE MM ENERGY
		
		puts "CALCULATING THE MM ENERGY FOR RESIDUE"

		set mm_energy_1 [mm_energy "pl9_only.pdb"]

		# CALCULATING THE QM ENERGY

		#puts "CALCULATING THE QM ENERGY FOR RESIDUE"
		
		#set qm_energy_1 [qm_energy "pl9_only.pdb"]
		
		set h1 [open "dummy" "r"]
		set data2 [read $h1]
		close $h1

		set k1 0
		set k3 0

		while { $k1 < [llength $data2] }	{
			set x [lindex $data2 $k1 0]
			set y [lindex $data2 $k1 1]
			set z [lindex $data2 $k1 2]
		
			set s [lindex $data1 $k1]
			puts "		**** FINDING THE $s NEAREST SOLVENT MOLECULE IN FRAME $fr ****"
			set count 0

			for {set j 0} {$j < $s} {incr j} {
				incr fd
				set k2 0
				set d 1000
				while {[lindex $data $k2] != "END" } {
					if { [lindex $data $k2] == "ATOM" && [lindex $data [expr { $k2 + 2 }]] == "[input13]" && [lindex $data [expr { $k2 + 3 }]] == "[input12]" } { 
						for {set px -$box_x} {$px <= $box_x} {set px [expr { $px + $box_x }]} { 
							for {set py -$box_y} {$py <= $box_y} {set py [expr { $py + $box_y }]} { 
								for {set pz -$box_z} {$pz <= $box_z} {set pz [expr { $pz + $box_z }]} { 

 									set coordx [lindex $data [expr { $k2 + 5 }]] 
									set coordx [expr { $coordx + $px }]
									set coordy [lindex $data [expr { $k2 + 6 }]] 
									set coordy [expr { $coordy + $py }]
									set coordz [lindex $data [expr { $k2 + 7 }]]
									set coordz [expr { $coordz + $pz }]
									set an1 [lindex $data [expr { $k2 + 4 }]]

									set delx [expr { $coordx - $x }]
									set dely [expr { $coordy - $y }]
									set delz [expr { $coordz - $z }]

									set delx
	
									set del($k2) [expr { sqrt(($delx*$delx) + ($dely*$dely) + ($delz*$delz)) }]

									if { $del($k2) < $d } {
										set count1 0
										set count2 0
										for {set l 0} {$l < $j} {incr l} {
											if { $k2 == $nn($l) } {
												incr count1
											}
										}
										for {set k4 0} {$k4 < $k3} {incr k4} {
											if { $rep($k4) == $an1  } {
												incr count2
											}
										}
										if { $count1 == 0 && $count2 == 0 } {
											set nn($j) $k2
											set d $del($k2)
											set xshift $px
											set yshift $py
											set zshift $pz
										}
									}
								}
							}
						}
					}
				incr k2
				}

				set k2 $nn($j)

				set next [expr { (1 - [input14])*11 }]

				for {set nsol 0} {$nsol < [input15]} {incr nsol} {

					set rn($nsol) [lindex $data [expr { $k2 + 1 + $next }]]
					set srn($nsol) [string length $rn($nsol)]

					set pdn($nsol) [lindex $data [expr { $k2 + 2 + $next }]]
					set spdn($nsol) [string length $pdn($nsol)]

					set an($nsol) [lindex $data [expr { $k2 + 4 + $next}]]
					set san($nsol) [string length $an($nsol)]

					set cordx($nsol) [lindex $data [expr { $k2 + 5 + $next }]]
					set cordx($nsol) [expr { $cordx($nsol) + $xshift }]
					set cordx($nsol) [format "%.3f" $cordx($nsol)]
					set scordx($nsol) [string length $cordx($nsol)]
		 
					set cordy($nsol) [lindex $data [expr { $k2 + 6 + $next}]] 
					set cordy($nsol) [expr { $cordy($nsol) + $yshift }]
					set cordy($nsol) [format "%.3f" $cordy($nsol)]
					set scordy($nsol) [string length $cordy($nsol)]

					set cordz($nsol) [lindex $data [expr { $k2 + 7 +$next}]] 
					set cordz($nsol) [expr { $cordz($nsol) + $zshift }]
					set cordz($nsol) [format "%.3f" $cordz($nsol)]
					set scordz($nsol) [string length $cordz($nsol)]

					set atom($nsol) [lindex $data [expr { $k2 + 10 + $next}]] 
					set satom($nsol) [string length $atom($nsol)]

					puts $h "ATOM  $p1($srn($nsol))$rn($nsol)$p($spdn($nsol))$pdn($nsol)  [input12] $p1($san($nsol))$an($nsol)     $c($scordx($nsol))$cordx($nsol) $c($scordy($nsol))$cordy($nsol) $c($scordz($nsol))$cordz($nsol)  1.00  0.00        $p($satom($nsol))$atom($nsol)" 

					incr next 11
				}
				puts $h "TER"	
				set rep($k3) $an(0)
				incr k3	
			}
		incr k1
		}
		close $h

		# CALCULATING THE MM ENERGY

		puts "CALCULATING THE MM ENERGY FOR SURROUNDING WATER"
		
		set mm_energy_2 [mm_energy "nwat.pdb"]
		puts "$mm_energy_2"

		# CALCULATING THE QM ENERGY

		#puts "CALCULATING THE QM ENERGY FOR SURROUNDING WATER"
		
		#set qm_energy_2 [qm_energy "nwat.pdb"]

		# COMBINING STRING 1 AND STRING 2

		set c1 [open "pl9_only.pdb" "r"]
		set file1 [read $c1]
		close $c1

		set c2 [open "nwat.pdb" "r"]
		set file2 [read $c2]
		close $c2

		#set file3 [concat $file1 $file2]
		set c3 [open "both.pdb" "w"]
		#puts $c3 "$file3"
		puts $c3 "$file1$file2"
		#puts $c3 "$file2"
		close $c3

		# CALCULATING THE MM ENERGY

		puts "CALCULATING THE MM ENERGY FOR RESIDUE + WATER"

		set mm_energy_3 [mm_energy "both.pdb"]

		# CALCULATING THE QM ENERGY

		puts "CALCULATING THE QM ENERGY FOR RESIDUE WATER AND RESIDUE + WATER USING COUNTERPOISE CALCULATION"

		qm_energy "both.pdb"

		set dum1 [open "dummy1" "r"]
		set dumm1 [read $dum1]
		close $dum1
		
		set qm_energy_1 [lindex $dumm1 0]
		set qm_energy_2 [lindex $dumm1 1]
		set qm_energy_3 [lindex $dumm1 2]

		set hy_en_mm [expr { $mm_energy_3 - $mm_energy_2 - $mm_energy_1 }]
		set mm_monomer [expr { $mm_energy_1 + $mm_energy_2 }]

		puts "$hy_en_mm		$qm_energy_3"
		#puts "$hy_en_mm"

		puts $mm "$fr	$mm_monomer	$mm_energy_3	$hy_en_mm"
	
		puts $qm "$fr	$qm_energy_1	$qm_energy_2	$qm_energy_3"
	}
	close $mm
	close $qm
}
	
proc mm_energy {a} {

	# FORMING THE STRIPPED PRMTOP FILE
	
	set inp [open "ai.dat" "w"]
	puts $inp "source leaprc.ff99SB"
	puts $inp "source leaprc.gaff"
	puts $inp "loadamberparams [input8]"
	if { [input16] > 0 } {
		puts $inp "loadamberparams [input17]"
		puts $inp "loadamberprep [input18]"
	}
	puts $inp "loadoff [input9]"
	puts $inp "pl9 = loadpdb $a"
	puts $inp "saveamberparm pl9 both.prmtop both.inpcrd"
	puts $inp "savepdb pl9 2.pdb"
	puts $inp "quit"
	close $inp

	set inp [open "ai.dat" "r"]

	exec tleap -f ai.dat

	close $inp

	# CALCULATING THE AMBER ENERGY FOR EACH FRAME

	# EXECUTING SANDER

	exec sander -O -i min.in -o min.out -p both.prmtop -c both.inpcrd -r min.rst

	set sander [open "min.out" "r"]
	set mm_en [read $sander]
	close $sander

	set var 0

	while { [lindex $mm_en $var] != "NSTEP" } {
		incr var
	}
	incr var 7

	set mm_energy [lindex $mm_en $var]


	return $mm_energy
}

proc qm_energy {a} {

	# GENERATING THE INPUT FILE USING ANTECHAMBER

	exec antechamber -i $a -fi pdb -o 1.com -fo gcrt -a both.prmtop -ao prmtop -ch frame.chk -gk "# hf/3-21G" -gm "%mem=10GB" -gn "%nprocshared=8"

	# FORMING THE INPUT FILE FOR COUNTERPOISE CALCULATION

	set h [open "1.com" "r"]
	set data [read $h]
	close $h

	set natoms [input2]

	set g [open "frame.com" "w"]

	puts $g "%nprocshared=4"
	puts $g "%mem=10GB"
	puts $g "%chk=frame.chk"
	puts $g "#p [input10] Counterpoise=2"
	puts $g ""
	puts $g "PER FRAME ENERGY CALCULATION"
	puts $g ""
	puts $g "0 1 0 1 0 1"

	set k 0

	while { [lindex $data $k] != 0 || [lindex $data [expr { $k + 1 }]] != 1 } {
		incr k
	}
	incr k 2

	set count $k
	set k1 $k
	while { $k < [llength $data] } {
		if { $count < [expr { $natoms + $k1 }] } {
			puts $g "	[lindex $data $k](Fragment=1)	[lindex $data [expr { $k + 1 }]]	[lindex $data [expr { $k + 2 }]]	[lindex $data [expr { $k + 3}]]"
			incr count
		} else { 
			puts $g "	[lindex $data $k](Fragment=2)	[lindex $data [expr { $k + 1 }]]	[lindex $data [expr { $k + 2 }]]	[lindex $data [expr { $k + 3}]]"
		}
		incr k 4
	}
	puts $g ""
	close $g

	exec g09 frame.com

	set f [open "frame.log" "r"]
	set data1 [read $f]
	close $f

	set dum [open "dummy1" "w"]

	set k 0
	
	while { [lindex $data1 $k] != "Counterpoise" || [lindex $data1 [expr { $k + 1 }]] != "corrected" || [lindex $data1 [expr { $k + 2 }]] != "energy"} {
		incr k
	}

	incr k 4

	set E2 [expr { [lindex $data1 $k] * 627.5095 }]

	incr k 9

	set E1 [expr { [lindex $data1 $k] * 627.5095 }]
	
	incr k 10

	set E3 [lindex $data1 $k] 
	
	puts $dum "$E1	$E2	$E3"
	close $dum

	#file delete frame.log

	after 1000
			
}

proc elem_cord {} {

	package require math::linearalgebra

	# GETTING THE COORDINATES WITH RESPECT TO CENTRE OF MASS EACH RESIDUE

	set f [open "frame_800.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "input1" "r"]
	set data1 [read $g]
	close $g

	# GEOMETRIX CENTER OF THE RING

	set gmx 0.0
	set gmy 0.0
	set gmz 0.0
	
	# HEAD POINT OF HEIGHT VECTOR

	set k1 0
	while { $k1 < [llength [lindex $data1 0]] } {
		set k2 0
		set atom [lindex $data1 0 $k1]
		while {[lindex $data $k2] != "TER" || [lindex $data [expr { $k2 - 1 }]] != "H" } {
			if { [lindex $data $k2] == "ATOM" &&  [lindex $data [expr { $k2 + 2 }]] == $atom } { 
				set coordx [lindex $data [expr { $k2 + 5 }]] 
				set coordy [lindex $data [expr { $k2 + 6 }]] 
				set coordz [lindex $data [expr { $k2 + 7 }]]
				set gmx [expr { $gmx + $coordx }]
				set gmy [expr { $gmy + $coordy }]
				set gmz [expr { $gmz + $coordz }]
			}
		incr k2
		}
	incr k1
	}
	set gmx [expr { $gmx / $k1 }]
	set gmy [expr { $gmy / $k1 }]
	set gmz [expr { $gmz / $k1 }] 
	puts "$gmx $gmy $gmz"

	# TAIL POINT OF HEIGHT VECTOR

	set gmx1 0.0
	set gmy1 0.0
	set gmz1 0.0

	set k1 0
	while { $k1 < [llength [lindex $data1 1]] } {
		set k2 0
		set atom [lindex $data1 1 $k1]
		while {[lindex $data $k2] != "TER" || [lindex $data [expr { $k2 - 1 }]] != "H" } {
			if { [lindex $data $k2] == "ATOM" &&  [lindex $data [expr { $k2 + 2 }]] == $atom } { 
				set coordx [lindex $data [expr { $k2 + 5 }]] 
				set coordy [lindex $data [expr { $k2 + 6 }]] 
				set coordz [lindex $data [expr { $k2 + 7 }]]
				set gmx1 [expr { $gmx1 + $coordx }]
				set gmy1 [expr { $gmy1 + $coordy }]
				set gmz1 [expr { $gmz1 + $coordz }]
			}
		incr k2
		}
	incr k1
	}
	set gmx1 [expr { $gmx1 / $k1 }]
	set gmy1 [expr { $gmy1 / $k1 }]
	set gmz1 [expr { $gmz1 / $k1 }] 
	puts "$gmx1 $gmy1 $gmz1"

	set gmx2 [expr { $gmx1 -$gmx }]
	set gmy2 [expr { $gmy1 -$gmy }]
	set gmz2 [expr { $gmz1 -$gmz }]

	set hvector [list [expr { $gmx1 - $gmx }] [expr { $gmy1 - $gmy }] [expr { $gmz1 - $gmz }]]

	set mod_hvector [expr { sqrt( ($gmx2*$gmx2) + ($gmy2*$gmy2) + ($gmz2*$gmz2) ) }]

	set unit_hvector [::math::linearalgebra::unitLengthVector $hvector]

	# GETTING THE WATER MOLECULE ON A CYLINDRICAL SHELL 10 A AWAY FROM THE MOLECULE

	set h [open "water_10A.pdb" "w"]

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""

	set c(4) "   "
	set c(5) "  "
	set c(6) " "
	set c(7) ""

	set k2 0
	set count 0	
	while {[lindex $data $k2] != "END" || [lindex $data [expr { $k2 - 1 }]] != "H" } {
		if { [lindex $data $k2] == "ATOM" && [lindex $data [expr { $k2 + 2 }]] == "O" && [lindex $data [expr { $k2 + 3 }]] == "WAT" } { 
			set coordx [lindex $data [expr { $k2 + 5 }]] 
			set coordy [lindex $data [expr { $k2 + 6 }]] 
			set coordz [lindex $data [expr { $k2 + 7 }]]
			set nx [expr { $coordx - $gmx }]
			set ny [expr { $coordy - $gmy }]
			set nz [expr { $coordz - $gmz }]
			set pvector [list $nx $ny $nz]
			set projection [::math::linearalgebra::dotproduct $hvector $pvector]
			set projection [expr { $projection / $mod_hvector }]
			set projection_vect [::math::linearalgebra::scale $projection $unit_hvector]
			set perpendicular_vect [::math::linearalgebra::sub $pvector $projection_vect]
			set perpendicular [::math::linearalgebra::norm $perpendicular_vect]

			set rn1 [lindex $data [expr { $k2 + 1 }]]
			set srn1 [string length $rn1]

			set an1 [lindex $data [expr { $k2 + 4 }]]
			set san1 [string length $an1]

			set rn2 [lindex $data [expr { $k2 + 12 }]]
			set srn2 [string length $rn2]

			set an2 [lindex $data [expr { $k2 + 15 }]]
			set san2 [string length $an2]

			set rn3 [lindex $data [expr { $k2 + 23 }]]
			set srn3 [string length $rn3]

			set an3 [lindex $data [expr { $k2 + 26 }]]
			set san3 [string length $an3]

			set x1 [lindex $data [expr { $k2 + 5 }]]
			set sx1 [string length $x1]
 
			set x2 [lindex $data [expr { $k2 + 16 }]]
			set sx2 [string length $x2]

			set x3 [lindex $data [expr { $k2 + 27 }]]
			set sx3 [string length $x3]

			set y1 [lindex $data [expr { $k2 + 6 }]] 
			set sy1 [string length $y1]

			set y2 [lindex $data [expr { $k2 + 17 }]]
			set sy2 [string length $y2]

			set y3 [lindex $data [expr { $k2 + 28 }]]
			set sy3 [string length $y3]

			set z1 [lindex $data [expr { $k2 + 7 }]] 
			set sz1 [string length $z1]

			set z2 [lindex $data [expr { $k2 + 18 }]]
			set sz2 [string length $z2]

			set z3 [lindex $data [expr { $k2 + 29 }]]
			set sz3 [string length $z3]


			if { $perpendicular > 10.0 && $perpendicular < 12.0 } {
				puts $h "ATOM   $p($srn1)$rn1  O   WAT  $p($san1)$an1     $c($sx1)$x1 $c($sy1)$y1 $c($sz1)$z1  1.00  0.00           O"  
				puts $h "ATOM   $p($srn2)$rn2  H1  WAT  $p($san2)$an2     $c($sx2)$x2 $c($sy2)$y2 $c($sz2)$z2  1.00  0.00           H"   
				puts $h "ATOM   $p($srn3)$rn2  H2  WAT  $p($san3)$an3     $c($sx3)$x3 $c($sy3)$y3 $c($sz3)$z3  1.00  0.00           H"    
				puts $h "TER"
				puts "$perpendicular"	
				incr count
			}		
		}
	incr k2
	}
puts "$count"
close $h
}

proc com {} {
	# GETTING THE COORDINATES WITH RESPECT TO CENTRE OF MASS EACH RESIDUE

	set f [open "1.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "parm" "r"]
	set data1 [read $g]
	close $g

	set k2 0

	set i 0
	set t 0

	set k1 0

	# READING THE MASS FROM THE PRMTOP FILE

	while { [lindex $data1 $k1] != "MASS" } {
		incr k1
	}
	incr k1 2

	# CENTRE OF MASS COORDINATES

	set comx 0.0
	set comy 0.0
	set comz 0.0
	set total_mass 0.0

	while {[lindex $data $k2] != "TER" || [lindex $data [expr { $k2 + 1 }]] != "ATOM" } {
		if { [lindex $data $k2] == "ATOM" } { 
			set mass [lindex $data1 $k1]
			set total_mass [expr {$total_mass + $mass}]
			incr k1
			set coordx [lindex $data [expr { $k2 + 5 }]] 
			set coordy [lindex $data [expr { $k2 + 6 }]] 
			set coordz [lindex $data [expr { $k2 + 7 }]]
			set comx [expr { $comx + ($mass*$coordx) }]
			set comy [expr { $comy + ($mass*$coordy) }]
			set comz [expr { $comz + ($mass*$coordz) }]
		}
		incr k2
	}
	set comx [expr { $comx / $total_mass }]
	set comy [expr { $comy / $total_mass }]
	set comz [expr { $comz / $total_mass }] 
	puts "$comx $comy $comz"
	return [list $comx $comy $comz]
}

proc diff_stages {} {

	package require math::linearalgebra

	set f [open "1.pdb" "r"]
	set data [read $f]
	close $f

	set centre [com]

	set g [open "dummy" "w"]
	
	puts $g "centre"

	set g [open "dummy" "r]
	set data1 [read $g]
	close $g

	set comx [lindex $data1 0]
	set comy [lindex $data1 1]
	set comz [lindex $data1 2]
	
	set k 0

	while { [lindex $data $k] != "END" || [lindex $data [expr {$k - 1}]] != "TER" } {
		if { [lindex $data $k] == "ATOM" && [lindex $data [expr { $k + 3 }]] == "WAT" } {
			set coordx [lindex $data [expr { $k + 5 }]]
			set coordx [lindex $data [expr { $k + 6 }]]
			set coordx [lindex $data [expr { $k + 7 }]]

			set delx [expr { $coordx - $comx }]
			set dely [expr { $coordy - $comy }]
			set delz [expr { $coordz - $comz }]

			set delr [list $delx $dely $delz]
		}
	}
}

# READING INPUT FILE

proc input1 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 0 1]]
}

proc input2 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 1 1]]
}

proc input3 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 2 1]]
}

proc input4 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 2 2]]
}

proc input5 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 2 3]]
}

proc input6 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 3 1]]
}

proc input7 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 4 1]]
}

proc input8 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 5 1]]
}

proc input9 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 5 2]]
}

proc input10 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 6 1]]
}

proc input11 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 7 1]]
}

proc input12 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 8 1]]
}

proc input13 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 8 2]]
}

proc input14 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 8 3]]
}
proc input15 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 8 4]]
}

proc input16 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	set a [llength [lindex $data 9]]

	return [expr { $a -1 }]
}

proc input17 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 9 1]]
}

proc input18 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 9 2]]
}

proc input19 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [list [lindex $data 1 2]]
}

#elem_cord
#com 
water_scale
nearest_wat


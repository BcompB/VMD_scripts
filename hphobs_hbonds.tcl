# Two procs are introduced here, hphobs to count hydrophobic interactions between two selections and hbonds to count H-bond interactions between two selections.
#  They can be used on the VMD command line or in tcl scripts after loading this.
# chmwvdk 2012

proc hphobs { seltext1 seltext2 name {dist 4.6} {timestep 1}} {

  puts "Usage: hphobs <\"selection 1\"> <\"selection 2\"> <name  for output file, i.e. hphobs_name.dat> <distance cutoff, default: 4.6> <timestep, default: 1>"
  puts "Running......"

# Selection text for heavy atoms that should be considered forming hydrophobic contacts
# Currently includes side chains of purely hydrophobic residues and Cs in THR,TYR,TRP that are NOT connected to O or N.
  set seltext_hphob "((hydrophobic and not resname TRP) or (resname TYR and not (name CZ or name OH)) or (resname THR and not (name CB or name OG1)) or (resname TRP and not (name NE1 or name CD1 or name CE2))) and not hydrogen and not (name N or name O or name C or name CA)"

  # start writing output file
  set outfile [open "hphobs_$name.dat" w]
  puts $outfile "#Number of hydrophobic interactions between \"$seltext1\" and \"$seltext2\" using $dist as distance cutoff"
  puts $outfile "#time(ps)\thphobs"

  # Loop over all frames
  set n [molinfo top get numframes]
  for { set i 0 } { $i < $n } { incr i } {
    set time [expr $i*$timestep]
    animate goto $i
    # initialize count
    set count 0
    # initialize hphob_list
    #??? set hphob_list [list] 
    # indexed list of all atoms involved in hphob interactions in selection 1
    set index_sel1 [[atomselect top "(within $dist of ($seltext2 and $seltext_hphob)) and ($seltext1 and $seltext_hphob)"] get index ]
    # indexed list of all atoms involved in hphob interactions in selection 2
    set index_sel2 [[atomselect top "(within $dist of ($seltext1 and $seltext_hphob)) and ($seltext2 and $seltext_hphob)"] get index ]

    # Loop over both lists, measure atom-atom distance and increase count if <= cutoff distance
    foreach at_sel1 $index_sel1 {
        foreach at_sel2 $index_sel2 {
            set checkdist [measure bond [list $at_sel1 $at_sel2]]
            if {$checkdist <= $dist} {incr count 1
            # insert code that gets atomname+resname+resid in contact and adds this to a list (to print in 'verbose' mode)
	    #lassign [$at_sel1 get {name resname resid}] name1 resname1 resid1
            #lassign [$at_sel2 get {name resname resid}] name2 resname2 resid2
            #lappend hphob_list [format "%s%s_%s---%s%s_%s" $resname1 $resid1 $name1 $resname2 $resid2 $name2] 
	    }
        }
    }

    # write count to output
    puts $outfile [format "%8d\t%d" $time $count]
  }

  close $outfile
  return
}


proc hbonds { seltext1 seltext2 name {dist 3.4} {angle 45} {timestep 1}} {
  puts "Usage: hbonds <\"selection 1\"> <\"selection 2\"> <name for output file, i.e. hbonds_name.dat> <distance cutoff, default: 3.4> <angle cutoff, default 45> <timestep, default: 1>"
  puts "Running......"

  set sel1 [atomselect top "($seltext1) and not name \"C.*\""]
  set sel2 [atomselect top "($seltext2) and not name \"C.*\""]

  set outfile [open "hbonds_$name.dat" w]
  puts $outfile "#Number of H-bonds between \"$seltext1\" and \"$seltext2\" using DA<=$dist, DHA +/- $angle"
  puts $outfile "#time(ps)\thbonds"
  set n [molinfo top get numframes]

  for { set i 0 } { $i < $n } { incr i } {
    set time [expr $i*$timestep]
    animate goto $i
    set hbonds_sel1 [llength [lindex [measure hbonds $dist $angle $sel1 $sel2] 0]]
    set hbonds_sel2 [llength [lindex [measure hbonds $dist $angle $sel2 $sel1] 0]]
    set hbonds [expr {$hbonds_sel1 + $hbonds_sel2}]
    puts $outfile [format "%8d\t%d" $time $hbonds]
  }

  close $outfile
  return
}


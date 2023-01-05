if [ "$1" = --description ] ; then
  echo "create a new collective variable from a template"
  exit 0
fi

if [ "$1" = --help ] || [ "$1" = -h ] ; then
  cat <<EOF
Usage:

  plumed partial_tempering scale_intra scale_inter < processed.top

where scale_intra and scale_inter are the Hamiltonian scaling factors and
processed.top is a post-processed topology file (i.e. produced with grompp -pp)
where each "hot" atom has a "_" appended to the atom type, e.g.:

     1 amber99_43_     1    RC5    O5'      1    -0.6223         16   ; qtot -0.6223

Notice that the section that should be edited is the [atoms] section for all the
molecules that you wish to affect (typically only for the solute, but you may also
want to change solvent parameters).

Also remember to first produce the processed.top file with grompp -pp. Editing a normal
topol.top file will not work, because it does not contain all the parameters.
The processed.top file should not have any "#include" statement.

# produce a processed topology
grompp -pp
# choose the "hot" atoms
vi processed.top
# generate the actual topology
plumed partial_tempering \$scale_intra \$scale_inter < processed.top > topol\$i.top

WARNING: It's not very robust and there might be force-field dependent issues!
A few tests are strongly suggested.

1. Compare partial_tempering with scale_intra=1.0 scale_inter=1.0 to non-scaled force field. E.g.
grompp -o topol-unscaled.tpr
grompp -pp
vi processed.top # choose the "hot" atoms appending "_". You can choose whatever.
plumed partial_tempering 1.0 1.0 < processed.top > topol-scaled.top # scale with factor 1
grompp -p topol-scaled.top -o topol-scaled.tpr
# Then do a rerun on a trajectory
mdrun -s topol-unscaled.tpr -rerun rerun.trr
mdrun -s topol-scaled.tpr -rerun rerun.trr
# and compare the resuling energy files. they should be identical

2. Compare partial_tempering with scale_intra=0.5 scale_inter=0.5 to non-scaled force field.
Repeat the same procedure but using "plumed partial_tempering 0.5 0.5".
Choose all the atoms in all the relevant [atoms] sections (e.g. solute, solvent and ions).
In the two resulting energy files you should see:
long range electrostatics, LJ, and dihedral energy is *half* in the scaled case
all other terms (bonds/bends) are identical.

EOF
  exit
fi

awk -v scale_intra=$1 -v scale_inter=$2 '
BEGIN{
  combrule=1;
  uniq=0;
}
function recname()
{
     if($1=="[" && $3=="]") return $2;
     return "";
}
function error(msg)
{
     print "ERROR:",msg > "/dev/stderr" ;
     exit;
}
function warning(msg)
{
     print "WARNING:",msg | "cat 1>&2"
}
{
# This is the suffix for "hot" atoms:
  suffix="_";

# format for writing parameters
  CONVFMT="%.12g"

##### PARSING DATABASE #####
# store comments:
  comments="";
# xrliu: match returns the character position at which that substring begins
# xrliu: substr returns substring starting at position a
  if(a=match($0,";")) comments=substr($0,a);
# remove comments:
  gsub(";.*","");
# echo empty line
  if(NF==0){
    print comments;
    next;
  }
# set name of current block
  if(recname() ) rec=recname();
# set defaults for nb interactions
  if(rec=="defaults" && NF==5) combrule=$2;
# if we are in atomtypes section, check which fields are present
# use same heuristics as in src/kernel/toppush.c
  if(rec=="atomtypes" && NF>=4){
    if((length($4)==1 && $4~"[a-zA-Z]")){
      bondtypefield=1;
      epsilonfield=6;
    }
    else if((length($6)==1 && $6~"[a-zA-Z]")){
      bondtypefield=2;
      epsilonfield=8;
    }
    else if((length($5)==1 && $5~"[a-zA-Z]")){
      if(substr($2,0,1) ~"[a-zA-Z]"){
        bondtypefield=2;
        epsilonfield=7;
      } else {
        # xrliu: this is the current format
        bondtypefield=1;
        sigmafield=6;
        epsilonfield=7;
      }
    } else error("in atomtypes");
    if(epsilonfield!=NF) error("in atomtypes");
# NOTE: OPLS uses bond types, thus we have to save the bondtype
# For other force fields (e.g. AMBER) atomtype is used as bondtype
# and column two is ignored (it is just atomic number).
    bondtype[$1]=$bondtypefield;
    sigma[$1]=$sigmafield;
    epsilon[$1]=$epsilonfield;
  }

# storing dihedraltypes:
  if(rec=="dihedraltypes" && ($5==1 || $5==2 || $5==3 || $5==4 || $5==5 || $5==9)){
    if($5==1 || $5==4 || $5==9) string=":"$6" "$7" "$8;
    else if($5==2) string=":"$6" "$7;
    else if($5==3) string=":"$6" "$7" "$8" "$9" "$10" "$11;
    else if($5==5) string=":"$6" "$7" "$8" "$9;
    type=$1"-"$2"-"$3"-"$4"-"$5
    if(type in params) params[type]=params[type]string; # xrliu: append string to params[type]
    else               params[type]=NR""string;
# parameters are commented since they are used inline
    print "; LINE("NR")",$0,comments
    # xrliu: The next statement forces awk to immediately stop processing the current record and go on to the next record
    next;
  }

##### SCANNING #####
# in case new list of atoms for a new molecule, delete the present list
  if(recname()=="atoms"){
    delete list_of_atoms;
    n_of_atoms=0;
  }
# detect amber type of each atom
  if(rec=="atoms" && NF>6){
     name=$2;
     gsub(suffix"$","",name)
     ato[$1]=name;
  }

##### PRINTING #####
# DIHEDRALS
  if(rec=="dihedrals" && ($5==1 || $5==2 || $5==3 || $5==4 || $5==5 || $5==9) ){
    found1=0; found4=0;
    for(j=0;j<n_of_atoms;j++) {
      if($1==list_of_atoms[j]) found1=1;
      if($4==list_of_atoms[j]) found4=1;
    }
    sscale=1.0;
    # xrliu: scale dihedral potentials using variable sscale
    if(found1)sscale*=sqrt(scale_intra);
    if(found4)sscale*=sqrt(scale_intra);

# this is the case in which dihedrals are online:
     if(NF>5){
       printf($1" "$2" "$3" "$4" "$5" ");
       if($5==1 || $5==4 || $5==9){
                                   if(NF!=8) error("dihedrals with type 1,4,9 should have 8 fields");
                                   printf($6" "$7*sscale" "$8);
       } else if($5==2) {
                                   if(NF!=7) error("dihedrals with type 2 should have 7 fields");
                                   printf($6" "$7*sscale);
       } else if($5==3) {
                                   if(NF!=11) error("dihedrals with type 3 should have 11 fields");
                                   printf($6*sscale" "$7*sscale" "$8*sscale" "$9*sscale" "$10*sscale" "$11*sscale);
       } else if($5==5) {
                                   if(NF!=9) error("dihedrals with type 5 should have 9 fields");
                                   printf($6*sscale" "$7*sscale" "$8*sscale" "$9*sscale);
       } else error("dihedrals with more than 5 fields should be 1,2,3,4,5 or 9");
       printf(" "comments"\n");
# this is the case in which we have to search the database
     } else if(NF==5){
       param="";
       atype[1]=bondtype[ato[$1]]
       atype[2]=bondtype[ato[$2]]
       atype[3]=bondtype[ato[$3]]
       atype[4]=bondtype[ato[$4]]
       
       progression=NR
       for(iswitch=0;iswitch<32;iswitch++){
         if(iswitch%2==0){
           a1=atype[1]; a2=atype[2]; a3=atype[3]; a4=atype[4];
         } else {
           a1=atype[4]; a2=atype[3]; a3=atype[2]; a4=atype[1];
         }
         if(int(iswitch/2)%2==1) a1="X";
         if(int(iswitch/4)%2==1) a2="X";
         if(int(iswitch/8)%2==1) a3="X";
         if(int(iswitch/16)%2==1) a4="X";
         test=a1"-"a2"-"a3"-"a4"-"$5;
         if(test in params){
           split(params[test],array,":");
           if(array[1]<progression){
             progression=array[1];
             param=params[test];
           }
         }
       }
       
    n=split(param,array,":");
    if(n<=1) error("params not found "$1" "$2" "$3" "$4" "$5" "atype[1]" "atype[2]" "atype[3]" "atype[4]);
    if($5!=9 && n!=2){
# in case of multiple dihedrals !=9, all parameters should be the same, otherwise I suspect there is some problem
      for(i=3;i<=n;i++){
        if((array[i]-array[2])**2>1e-20) error("multiple dihedrals !=9: parameters "array[i]" and "array[2]" are different\n");
      }
# then, I just take one of the instances
      param=array[1]":"array[2];
      n=split(param,array,":");
    }
    
    printf("; parameters for types %s %s %s %s at LINE(%s)\n",atype[1],atype[2],atype[3],atype[4],array[1]);
    for(i=2;i<=n;i++){
      printf($1" "$2" "$3" "$4" "$5" ");
      split(array[i],array1," ");
      if($5==1 || $5==4 || $5==9){
                                  printf(array1[1]" "array1[2]*sscale" "array1[3]);
      } else if($5==2) {
                                  printf(array1[1]" "array1[2]*sscale);
      } else if($5==3) {
                                  printf(array1[1]*sscale" "array1[2]*sscale" "array1[3]*sscale" "array1[4]*sscale" "array1[5]*sscale" "array1[6]*sscale);
      } else if($5==5) {
                                  printf(array1[1]*sscale" "array1[2]*sscale" "array1[3]*sscale" "array1[4]*sscale);
      } else error("dihedrals with more than 5 fields should be 1,2,3,4,5 or 9");
      printf(comments);
      printf("\n");
    }
   } else error("dihedrals should have at least 5 fields");
# ATOMTYPES
  } else if(rec=="atomtypes" && NF>=4){
    for(i=1;i<NF-1;i++)printf($i" "); 
    if (combrule==2) {
       printf("%.7E %.7E %s \n",4*$NF*$(NF-1)**6,4*$NF*$(NF-1)**12,comments);
    } else if (combrule==1) {
       printf("%.7E %.7E %s \n", $(NF-1), $NF, comments);
    }
    # scale LJ
    printf($1""suffix" "bondtype[$1]" ");
    from=3;
    if(NF==6) from=2; # GROMOS does not store bondtype by default, so we should add one column
    for(i=from;i<NF-1;i++)printf($i" "); 
    if (combrule==2) {
       printf("%.7E %.7E %s \n",scale_intra*4*$NF*$(NF-1)**6,scale_intra*4*$NF*$(NF-1)**12, "; scaled");
    } else if (combrule==1) {
       printf("%.7E %.7E %s \n", scale_intra*$(NF-1),scale_intra*$NF, " ; scaled");
    }
# ATOMTYPES (PAIRS)
  } else if((rec=="pairtypes" || rec=="nonbond_params") && NF>=5){
    printf("%s %s %i %.7E %.7E ",$1,$2,$3,4*$5*$4**6,4*$5*$4**12)
    print comments;
    if (combrule==2) {
       printf("%s %s %i %.7E %.7E %s\n", $1""suffix,$2,$3,scale_inter*4*$5*$4**6,scale_inter*4*$5*$4**12," ; scaled");
       printf("%s %s %i %.7E %.7E %s\n", $1,$2""suffix,$3,scale_inter*4*$5*$4**6,scale_inter*4*$5*$4**12," ; scaled");
       printf("%s %s %i %.7E %.7E %s\n", $1""suffix,$2""suffix,$3,scale_intra*4*$5*$4**6,scale_intra*4*$5*$4**12," ; scaled");
    } else if (combrule==1) {
       printf("%s %s %i %.7E %.7E %s\n", $1""suffix,$2,$3,scale_inter*$4,scale_inter*$5," ; scaled");
       printf("%s %s %i %.7E %.7E %s\n", $1,$2""suffix,$3,scale_inter*$4,scale_inter*$5," ; scaled");
       printf("%s %s %i %.7E %.7E %s\n", $1""suffix,$2""suffix,$3,scale_intra*$4,scale_intra*$5," ; scaled");
    }
# ATOMS
  } else if(rec=="atoms" && NF>=7){
     if($2~".*"suffix"$"){
       # xrliu: scale solute charge
       if(NF>=8) print $1,$2,$3,$4,$5,$6,$7*sqrt(scale_intra),$8,comments;
       if(NF==7) print $1,$2,$3,$4,$5,$6,$7*sqrt(scale_intra),comments;
       list_of_atoms[n_of_atoms]=$1;
       # xrliu
       if ($2 in solute_atom_type) {
          solute_atom_type[$2]=solute_atom_type[$2]"-"$1;
       } else {
          solute_atom_type[$2]=$1;
          uniq++;
       }
       #print "xrliu debug:", uniq, $2, solute_atom_type[$2] > "/dev/stderr";
       n_of_atoms++;
     }
     else {
       print $0;
       if ($2 in solvent_atom_type) {
          solvent_atom_type[$2]=solvent_atom_type[$2]"-"$1;
       } else {
          solvent_atom_type[$2]=$1;
       }
     }
  } else if (rec=="defaults" && NF==5) {
    print("1 1",$3,$4,$5); 
# EVERYTHING ELSE (just print)
  } else print $0,comments
}
END {
print("");
print(";[ nonbond_params ] ;for solute-solvent interactions");
print("; i    j    funct    c6   c12");
for (typei in solute_atom_type) {
     typeii=substr(typei, 0, length(typei)-1);
     for (typej in solvent_atom_type) {
         sij=0.5*(sigma[typeii]+sigma[typej]);
         eij=sqrt(epsilon[typeii]*epsilon[typej]);
         printf("%s %s %i %.7E %.7E \n", typei, typej, 1, scale_inter*4*eij*sij**6, scale_inter*4*eij*sij**12);
     }
}
print(";[ nonbond_params ] ;for solute-solute interactions");
print(";number of solute atom types:",length(solute_atom_type));
i=0;
for (typei in solute_atom_type) {
     typeii=substr(typei, 0, length(typei)-1);
     j=0;
     for (typej in solute_atom_type) {
         typejj=substr(typej, 0, length(typej)-1);
         sij=0.5*(sigma[typeii]+sigma[typejj]);
         eij=sqrt(epsilon[typeii]*epsilon[typejj]);
         if (j>i) {
            printf("%s %s %i %.7E %.7E \n", typei, typej, 1, scale_intra*4*eij*sij**6, scale_intra*4*eij*sij**12);
            #printf("%s %s %i %.7E %.7E ", typei, typej, 1, scale_intra*4*eij*sij**6, scale_intra*4*eij*sij**12);
            #print(i,j);
         }
         j++;
     }
     i++;
}
}
'


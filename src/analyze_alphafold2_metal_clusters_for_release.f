
C this is code for release to the community
C for internal purposes note that it should be functionally identical
C to:

C analyze_alphafold2_metal_regions_final_check_clashes_general_extended_clashfix.f



      program fit_ligand_molecules
 
      implicit real(a-h,o-z)
      parameter   (pi=3.141592654)
      INTEGER, PARAMETER :: K4B=selected_int_kind(9)

      interface
      subroutine pdbsup(nmov_loc,nmtc_loc,xm_loc,ym_loc,zm_loc,
     &                  xf_loc,yf_loc,zf_loc,rms_loc)

      integer, intent(in)     :: nmov_loc
      integer, intent(in)     :: nmtc_loc
      real, intent(inout)     :: xm_loc(1:nmov_loc)
      real, intent(inout)     :: ym_loc(1:nmov_loc)
      real, intent(inout)     :: zm_loc(1:nmov_loc)
      real, intent(inout)     :: xf_loc(1:nmtc_loc)
      real, intent(inout)     :: yf_loc(1:nmtc_loc)
      real, intent(inout)     :: zf_loc(1:nmtc_loc)
      real,    intent(out)    :: rms_loc

      end subroutine
      end interface


      type gridtype
        integer                     :: num
        real, dimension(:), pointer :: x
        real, dimension(:), pointer :: y
        real, dimension(:), pointer :: z
        integer, dimension(:), pointer :: ires
        character*4, dimension(:), pointer :: anam
        character*3, dimension(:), pointer :: rnam
      end type gridtype

      type (gridtype), allocatable  :: grid(:,:,:)
     

C up to 1000 chains
C up to 10000 residues in each chain

      integer              :: idone(10000)
      integer              :: rnum_tmp
      character*4          :: anam_tmp
      character*3          :: rnam_tmp
 
      integer              :: rnum(10000)
      integer              :: natm(10000)
      character*4          :: anam(1000,10000)
      character*4          :: anam2(1000)
      character*3          :: rnam(10000)
      character*80         :: char80(1000,10000)
 
      character*100        :: initial_pdb

      character*4          :: char4

      real                 :: beta(100000)
      real xf(1000000) ! amino acids already there
      real yf(1000000)
      real zf(1000000)
      real xm(10000) ! amino acid to be added
      real ym(10000)
      real zm(10000)

C ligand stuff

      character*40         :: outfile1
      character*40         :: outfile2
      character*23         :: outfile3
      character*30         :: outfile4
      character*27         :: outfile5
      character*30         :: outfile6
      character*32         :: outfile7
      character*30         :: outfile8

      character*80 junk        
      character*80 string

      integer              :: ncl_disulf(1000)
      integer              :: idr_disulf(1000)
      integer              :: jdr_disulf(1000)
      real                 :: dis_disulf(1000)

      real                 :: rms_best(1000)
      real                 :: rms_typ(1000)
      integer              :: i_best(1000)
      integer              :: j_best(1000)
      integer              :: k_best(1000)
      integer              :: l_best(1000)
      integer              :: m_best(1000)
      integer              :: n_best(1000)
      character*6          :: char6_arr(1000)
      character*6          :: char6_brr(1000)
      character*6          :: char6 
      character*6          :: char6_job
      character*30         :: char30
      character*12         :: char12
      character*30         :: string1(1000,1000)
      character*12         :: string2(1000,1000)
      real                 :: rms_thresh(1000)
      integer              :: nmtc(1000)
      integer              :: nmov(1000)
      integer              :: ncys_lig(1000)
      integer              :: nhis_lig(1000)
      integer              :: nasp_lig(1000)
      integer              :: ndoe_lig(1000)
      integer              :: nbbn_lig(1000)
      real                 :: xlig(1000,1000)
      real                 :: ylig(1000,1000)
      real                 :: zlig(1000,1000)
      integer              :: ifound(1000) ! yes/no
      integer              :: nfound(1000) ! how many?
 
      real                 :: xfin(100000) 
      real                 :: yfin(100000) 
      real                 :: zfin(100000) 

      integer              :: my_unit1(20,20)
      integer              :: ispecial_region(10000)
      integer              :: num_mem_region(10000)
      integer              :: num_mem_region_orig(10000)
      integer              :: dun_mem_region(10000) ! =#ligsite if done
      integer              :: idr_mem_region(10000,1000) ! res#
      integer              :: nlg_mem_region(10000,1000) ! #poss ligs
      integer              :: mlg_mem_region(1000) ! #actl ligs
      integer              :: ilg_mem_region(1000,3) ! ligsite
      character*3          :: nam_mem_region(10000,1000) ! res name
      character*4          :: atm_mem_region(10000,1000) ! atm name
      integer              :: ida_mem_region(10000,1000) ! atm#
      real                 :: xxx_mem_region(10000,1000) ! xcoord
      real                 :: yyy_mem_region(10000,1000) ! ycoord
      real                 :: zzz_mem_region(10000,1000) ! zcoord
      real                 :: bfc_mem_region(10000,1000) ! bfactor

      integer              :: ires_tot_atm(100000)
      character*3          :: rnam_tot_atm(100000)
      character*4          :: anam_tot_atm(100000)
      real                 :: xxx_tot_atm(100000)
      real                 :: yyy_tot_atm(100000)
      real                 :: zzz_tot_atm(100000)

      integer              :: idr_mem_final(1000)
      character*3          :: nam_mem_final(1000)
      integer              :: idr_mem_tmp(1000)
      character*3          :: nam_mem_tmp(1000)

      character*80 ligand_pdb(1000)
      character*80 ligand_tmp
      character*80 ligand_list
      character*80 file_pdb

C read all command line arguments

      call getarg(1,junk)
      read(junk,*)file_pdb     ! input pdb file
      call getarg(2,junk)
      read(junk,*)ligand_list ! file listing ligands to try fitting
      call getarg(3,junk)
      read(junk,*)dist_lo  ! include in region if dist>=this
      call getarg(4,junk)
      read(junk,*)dist_hi  ! include in region if dist<=this
      call getarg(5,junk)
      read(junk,*)clash_dist ! steric clashes if dist <this
      call getarg(6,junk)
      read(junk,*)disulf_dist ! disulfide if dist<=this
      call getarg(7,junk)
      read(junk,*)bbnadd_dist ! add bbn atom to region if dist<=this
      call getarg(8,junk)
      read(junk,*)clash_lig_dist ! lig-lig clash dist
      call getarg(9,junk)
      read(junk,*)job_number ! use this to distinguish jobs

      disulf_dist2=disulf_dist**2
      bbnadd_dist2=bbnadd_dist**2
      clash_lig_dist2=clash_lig_dist**2

      write(char6_job,'(i6)')job_number
      do i=1,6
        if(char6_job(i:i).eq.' ') char6_job(i:i)='0'
      enddo

C outfile1 records, for each searched ligand, how many places we
C successfully placed it

C outfile2 records, for each searched ligand, *whether* we were able to
C successfully place it

C outfile3 records the superimposed coordinates of all placed ligands

C outfile4 records text information regarding every placed ligand

C outfile5 records all regions and members found in this pdb file -
C many of these will remain unfitted by ligands

C outfile6 records all putative disulfides in the protein

C outfile7 records, for each 4-CYS-only region, the RMSDs of all
C ligands that contain only 4CYSs

C outfile8 records all single-region CYS, all putative disulfides, and
C all orfan CYS residues that remain unassigned in regions

      outfile1='ligands_placed_number_summary_XXXXXX.txt'
      outfile2='ligands_placed_yes_no_summary_XXXXXX.txt'
      outfile3='ligand_sites_XXXXXX.pdb'
      outfile4='ligand_summary_info_XXXXXX.txt'
      outfile5='region_info_all_XXXXXX.txt'
      outfile6='putative_disulfides_XXXXXX.txt'
      outfile7='four_cys_ligand_RMSDs_XXXXXX.txt'
      outfile8='cys_singles_doubles_XXXXXX.txt'

      outfile1(31:36)=char6_job
      outfile2(31:36)=char6_job
      outfile3(14:19)=char6_job
      outfile4(21:26)=char6_job
      outfile5(17:22)=char6_job
      outfile6(21:26)=char6_job
      outfile7(23:28)=char6_job
      outfile8(21:26)=char6_job

      clash2=clash_dist**2

      write(*,*)
      write(*,*)file_pdb
      write(*,*)

      ires_mode=3 ! hardwire code to do CYS + HIS

      nfin_atm=0  ! #atoms to write to ligand_sites.pdb
 
      dist_lo2=dist_lo**2
      dist_hi2=dist_hi**2

C set to zero some later arrays

      dun_mem_region=0
      mlg_mem_region=0
      ilg_mem_region=0

C open up the file listing all ligand types

      nmtc_min=999999
      nmtc_max=-1
      nlig=0

C also need to check if we're looking for backbone N or O - currently we
C allow only *ONE* per run of the code - may need to extend this later

      iwantbbnN=0
      iwantbbnO=0

      open(unit=12,file=ligand_list,status='unknown')
501   read(12,*,end=502)ligand_tmp,ntmp,rtmp

      nlig=nlig+1

      ligand_pdb(nlig)=ligand_tmp

      nmtc(nlig)=ntmp ! set #atoms used in rmsd calculation here
      rms_thresh(nlig)=rtmp

      nmtc_min=min(nmtc_min,nmtc(nlig)) 
      nmtc_max=max(nmtc_max,nmtc(nlig)) 

C read in the ligand pdb file to identify all atoms of ligand
C the first nmtc atoms are the atoms used for the fitting and
C the remaining atoms are the actual region to be placed

C first figure out how many atom

        open(unit=11,file=ligand_pdb(nlig),status='unknown')
        ntmp=0
801     read(11,'(a4)',end=802)char4
        if(char4.eq.'ATOM') ntmp=ntmp+1
        goto 801
802     close(11)
        nmov(nlig)=ntmp
      
C now read all the coords etc, as well as the # of CYS and HIS residues
C that are to be used in the fit - if, e.g. we have 2CYS and 2HIS then
C we'll only try to fit it to region-atom combinations that contain
C 2CYS and 2HIS etc
C GENERAL - hardwire looking for others too

        ncys_lig(nlig)=0
        nhis_lig(nlig)=0
        nasp_lig(nlig)=0 ! ASP 
        ndoe_lig(nlig)=0 ! ASP or GLU
        nbbn_lig(nlig)=0 ! Backbone

        open(unit=11,file=ligand_pdb(nlig),status='unknown')
        do m=1,nmov(nlig)
          read(11,'(a30,3f8.3,a12)')char30,x,y,z,char12
          string1(nlig,m)=char30
          string2(nlig,m)=char12
          xlig(nlig,m)=x
          ylig(nlig,m)=y
          zlig(nlig,m)=z

C if this is a match atom also record the # of CYS and HIS

          if(m.le.nmtc(nlig)) then
            if(char30(18:20).eq.'CYS') ncys_lig(nlig)=ncys_lig(nlig)+1
            if(char30(18:20).eq.'HIS') nhis_lig(nlig)=nhis_lig(nlig)+1
            if(char30(18:20).eq.'ASP') nasp_lig(nlig)=nasp_lig(nlig)+1
            if(char30(18:20).eq.'DOE') ndoe_lig(nlig)=ndoe_lig(nlig)+1
            if(char30(18:20).eq.'BBN') then
              nbbn_lig(nlig)=nbbn_lig(nlig)+1
              if(char30(13:15).eq.' N ') iwantbbnN=1 
              if(char30(13:15).eq.' O ') iwantbbnO=1 
            endif
          endif

        enddo
        close(11)

C debug statement here

c       write(*,*)'LIGAND INFO ',nlig,nmtc(nlig),nmov(nlig),ligand_pdb,
c    &             ncys_lig(nlig),nhis_lig(nlig)

      goto 501
502   close(12)

C put in a check that we're not looking for *both* backbone N and Os

      if(iwantbbnN.eq.1.and.iwantbbnO.eq.1) then
        write(*,*)
        write(*,*)'FATAL ERROR - can only search BBN N or O'
        write(*,*)'cannot search for *both* N and O in same run'
        write(*,*)
        stop
      endif
        
C put in a check on nmtc_min and nmtc_max

      if(nmtc_min.lt.3) then
        write(*,*)
        write(*,*)'FATAL ERROR - #match atoms must be >=3'
        write(*,*)
        stop
      endif

      if(nmtc_max.gt.6) then
        write(*,*)
        write(*,*)'FATAL ERROR - #match atoms must be <=6'
        write(*,*)
        stop
      endif

C zero out ifound and nfound for each ligand type

      do n=1,nlig
        ifound(n)=0
        nfound(n)=0
      enddo

C now think about opening up the protein pdb file for analysis

C first let's figure out which residue types we're looking for:
C note that we don't look for BBN atoms just yet - we look for them
C *after* we've finalized our regions

      iwantcys=0
      iwanthis=0
      iwantasp=0
      iwantdoe=0
      iwantbbn=0

      do n=1,nlig
        if(ncys_lig(n).gt.0) iwantcys=1
        if(nhis_lig(n).gt.0) iwanthis=1
        if(nasp_lig(n).gt.0) iwantasp=1
        if(ndoe_lig(n).gt.0) iwantdoe=1
        if(nbbn_lig(n).gt.0) iwantbbn=1
      enddo

C set the total # of identified ligand binding sites to zero...

      nligsites=0
      nfails=0 ! # of non-assigned regions

C read in the pdb file to identify all the residues

      do n=1,30
        initial_pdb(n:n)=' '
      enddo

      len1=len(trim(file_pdb)) ! WTF this is overkill
      write(*,*)'working on ',file_pdb(1:len1)
ccc   write(*,*)'working on ',initial_pdb(1:len1)

      initial_pdb(1:len1)=file_pdb(1:len1)

      rnum_last=-999
      nres=0

C also store all the atom info for putting them on a grid

      num_tot_atm=0

C open up a file that stores all located ligands

      open(unit=31,file=outfile4,status="unknown")

C open up a file that stores all 4-CYS RMSDs

      open(unit=33,file=outfile7,status="unknown")

C set the number of regions to zero then add to them as we find new
C residues that match the resnames we're looking for - after we're done
C reading the pdb file we'll have an initial list of "regions" that
C we'll then seek to merge to form true regions

      num_region=0

      open(unit=11,file=initial_pdb(1:len1),status="unknown")
10    read(11,'(a80)',end=20)string

C skip non-ATOM entries

      if(string(1:4).ne.'ATOM') goto 10

C now read the entirety of the ATOM entry

      read(string,'(12x,a4,1x,a3,2x,i4,4x,3f8.3,2f6.2)')
     &     anam_tmp,rnam_tmp,rnum_tmp,xt,yt,zt,occ,beta_tmp

C store the atom infor for the grid

      num_tot_atm=num_tot_atm+1
      xxx_tot_atm(num_tot_atm)=xt
      yyy_tot_atm(num_tot_atm)=yt
      zzz_tot_atm(num_tot_atm)=zt
      ires_tot_atm(num_tot_atm)=rnum_tmp
      anam_tot_atm(num_tot_atm)=anam_tmp
      rnam_tot_atm(num_tot_atm)=rnam_tmp

C if this is a new residue number then increment nres etc

      if(rnum_tmp.ne.rnum_last) then

        nres=nres+1
        natm(nres)=1
        rnam(nres)=rnam_tmp
        rnum(nres)=rnum_tmp
        anam(natm(nres),nres)=anam_tmp
        char80(natm(nres),nres)=string
        beta(nres)=beta_tmp
        rnum_last=rnum_tmp

        if(rnam_tmp.eq.'CYS') num_cys_total=num_cys_total+1

C since this is the first atom of a new residue we can ask if this is
C one of the residues that we're following - if yes then increment the
C number of regions and set num_mem_region for this one to zero... we'll
C add to this when we read subsequent atoms of the same residue later

        if(iwantcys.eq.1) then
          if(rnam_tmp.eq.'CYS') then
            num_region=num_region+1
            num_mem_region(num_region)=0
          endif
        endif
        if(iwanthis.eq.1) then
          if(rnam_tmp.eq.'HIS') then
            num_region=num_region+1
            num_mem_region(num_region)=0
          endif
        endif
        if(iwantasp.eq.1) then
          if(rnam_tmp.eq.'ASP') then
            num_region=num_region+1
            num_mem_region(num_region)=0
          endif
        endif
        if(iwantdoe.eq.1) then
          if(rnam_tmp.eq.'ASP'.or.rnam_tmp.eq.'GLU') then
            num_region=num_region+1
            num_mem_region(num_region)=0
          endif
        endif

C if this is *not* the first atom of the residue (which, for alphafold2
C always appears to be a N) then just add this one to the residue 

      else

        natm(nres)=natm(nres)+1
        anam(natm(nres),nres)=anam_tmp
        char80(natm(nres),nres)=string

      endif

C now we also want to check if this is an atom to add to a region

        if(iwantcys.eq.1) then
          if(rnam_tmp.eq.'CYS') then
          if(anam_tmp(2:3).eq.'SG') then
          num_mem_region(num_region)=num_mem_region(num_region)+1
          idr_mem_region(num_region,num_mem_region(num_region))=nres
          nam_mem_region(num_region,num_mem_region(num_region))='CYS'
          atm_mem_region(num_region,num_mem_region(num_region))=anam_tmp
          xxx_mem_region(num_region,num_mem_region(num_region))=xt
          yyy_mem_region(num_region,num_mem_region(num_region))=yt
          zzz_mem_region(num_region,num_mem_region(num_region))=zt
          bfc_mem_region(num_region,num_mem_region(num_region))=beta_tmp

C set #possible ligands for this site (nlg)

          nlg_mem_region(num_region,num_mem_region(num_region))=2

          endif
          endif
        endif
        if(iwanthis.eq.1) then
          if(rnam_tmp.eq.'HIS') then
          if(anam_tmp(2:4).eq.'ND1'.or.
     &       anam_tmp(2:4).eq.'NE2') then
          num_mem_region(num_region)=num_mem_region(num_region)+1
          idr_mem_region(num_region,num_mem_region(num_region))=nres
          nam_mem_region(num_region,num_mem_region(num_region))='HIS'
          atm_mem_region(num_region,num_mem_region(num_region))=anam_tmp
          xxx_mem_region(num_region,num_mem_region(num_region))=xt
          yyy_mem_region(num_region,num_mem_region(num_region))=yt
          zzz_mem_region(num_region,num_mem_region(num_region))=zt
          bfc_mem_region(num_region,num_mem_region(num_region))=beta_tmp

C set #possible ligands for this site (nlg)

          nlg_mem_region(num_region,num_mem_region(num_region))=1

          endif
          endif
        endif
        if(iwantasp.eq.1) then
          if(rnam_tmp.eq.'ASP') then
          if(anam_tmp(2:4).eq.'OD1'.or.
     &       anam_tmp(2:4).eq.'OD2') then
          num_mem_region(num_region)=num_mem_region(num_region)+1
          idr_mem_region(num_region,num_mem_region(num_region))=nres
          nam_mem_region(num_region,num_mem_region(num_region))='ASP'
          atm_mem_region(num_region,num_mem_region(num_region))=anam_tmp
          xxx_mem_region(num_region,num_mem_region(num_region))=xt
          yyy_mem_region(num_region,num_mem_region(num_region))=yt
          zzz_mem_region(num_region,num_mem_region(num_region))=zt
          bfc_mem_region(num_region,num_mem_region(num_region))=beta_tmp
          nlg_mem_region(num_region,num_mem_region(num_region))=2
          endif
          endif
        endif
        if(iwantdoe.eq.1) then
          if(rnam_tmp.eq.'ASP') then
          if(anam_tmp(2:4).eq.'OD1'.or.
     &       anam_tmp(2:4).eq.'OD2') then
          num_mem_region(num_region)=num_mem_region(num_region)+1
          idr_mem_region(num_region,num_mem_region(num_region))=nres
          nam_mem_region(num_region,num_mem_region(num_region))='DOE'
          atm_mem_region(num_region,num_mem_region(num_region))=anam_tmp
          xxx_mem_region(num_region,num_mem_region(num_region))=xt
          yyy_mem_region(num_region,num_mem_region(num_region))=yt
          zzz_mem_region(num_region,num_mem_region(num_region))=zt
          bfc_mem_region(num_region,num_mem_region(num_region))=beta_tmp
          nlg_mem_region(num_region,num_mem_region(num_region))=2
          endif
          elseif(rnam_tmp.eq.'GLU') then
          if(anam_tmp(2:4).eq.'OE1'.or.
     &       anam_tmp(2:4).eq.'OE2') then
          num_mem_region(num_region)=num_mem_region(num_region)+1
          idr_mem_region(num_region,num_mem_region(num_region))=nres
          nam_mem_region(num_region,num_mem_region(num_region))='DOE'
          atm_mem_region(num_region,num_mem_region(num_region))=anam_tmp
          xxx_mem_region(num_region,num_mem_region(num_region))=xt
          yyy_mem_region(num_region,num_mem_region(num_region))=yt
          zzz_mem_region(num_region,num_mem_region(num_region))=zt
          bfc_mem_region(num_region,num_mem_region(num_region))=beta_tmp
          nlg_mem_region(num_region,num_mem_region(num_region))=2
          endif
          endif
        endif

      goto 10
20    close(11)

C now we're done reading the pdb file so we can write out our initial
C "region" definitions, i.e. all residues that we'll be tracking

      write(*,*)'just read #res = ',nres
      do n=1,num_region
        do m=1,num_mem_region(n)
          write(*,838)n,m,nam_mem_region(n,m),atm_mem_region(n,m),
     &                    idr_mem_region(n,m),nlg_mem_region(n,m),
     &                    xxx_mem_region(n,m),zzz_mem_region(n,m),
     &                    zzz_mem_region(n,m),bfc_mem_region(n,m)
        enddo
      enddo
838   format('initial region check ',2i6,1x,a3,1x,a4,1x,2i6,4f10.3)

C ----------------------------------------------------------------------
C now we need to assign all atoms to grid cells...

      xmin= 9999.9
      ymin= 9999.9
      zmin= 9999.9
      xmax=-9999.9
      ymax=-9999.9
      zmax=-9999.9

      do i=1,num_tot_atm
        xmin=min(xmin,xxx_tot_atm(i))
        ymin=min(ymin,yyy_tot_atm(i))
        zmin=min(zmin,zzz_tot_atm(i))
        xmax=max(xmax,xxx_tot_atm(i))
        ymax=max(ymax,yyy_tot_atm(i))
        zmax=max(zmax,zzz_tot_atm(i))
      enddo

      rcll=5.0 ! use a 5A grid
      xlen=xmax-xmin
      ylen=ymax-ymin
      zlen=zmax-zmin
      xinv=1.0/xlen
      yinv=1.0/ylen
      zinv=1.0/zlen
      numx=int(xlen/rcll)+1
      numy=int(ylen/rcll)+1
      numz=int(zlen/rcll)+1

      write(*,*)'atoms grid sizes ',numx,numy,numz
      write(*,*)'xmin xmax xlen ',xmin,xmax,xlen
      write(*,*)'ymin ymax ylen ',ymin,ymax,ylen
      write(*,*)'zmin zmax zlen ',zmin,zmax,zlen

      allocate(grid(1:numx,1:numy,1:numz))

      grid%num=0
      num_in_box_max=0

      do n=1,num_tot_atm

        i=int((xxx_tot_atm(n)-xmin)/rcll)+1
        j=int((yyy_tot_atm(n)-ymin)/rcll)+1
        k=int((zzz_tot_atm(n)-zmin)/rcll)+1

        if(i.le.0) i=numx
        if(j.le.0) j=numy
        if(k.le.0) k=numz
        if(i.gt.numx) i=1
        if(j.gt.numy) j=1
        if(k.gt.numz) k=1

        grid(i,j,k)%num=grid(i,j,k)%num+1
        if(grid(i,j,k)%num.gt.num_in_box_max)
     &     num_in_box_max=grid(i,j,k)%num

      enddo

      write(*,120)num_atm,num_in_box_max
120   format('just assigned ',i8,' atoms to the grid ',
     &       'max number in any one grid box = ',i8)

C now allocate the memory and fill them in

      do k=1,numz
        do j=1,numy
          do i=1,numx
            allocate(grid(i,j,k)%x(1:grid(i,j,k)%num))
            allocate(grid(i,j,k)%y(1:grid(i,j,k)%num))
            allocate(grid(i,j,k)%z(1:grid(i,j,k)%num))
            allocate(grid(i,j,k)%ires(1:grid(i,j,k)%num))
            allocate(grid(i,j,k)%anam(1:grid(i,j,k)%num))
            allocate(grid(i,j,k)%rnam(1:grid(i,j,k)%num))
            grid(i,j,k)%num=0
          enddo
        enddo
      enddo

      do n=1,num_tot_atm

        i=int((xxx_tot_atm(n)-xmin)/rcll)+1
        j=int((yyy_tot_atm(n)-ymin)/rcll)+1
        k=int((zzz_tot_atm(n)-zmin)/rcll)+1

        if(i.le.0) i=numx
        if(j.le.0) j=numy
        if(k.le.0) k=numz
        if(i.gt.numx) i=1
        if(j.gt.numy) j=1
        if(k.gt.numz) k=1

        grid(i,j,k)%num=grid(i,j,k)%num+1
        grid(i,j,k)%x(grid(i,j,k)%num)=xxx_tot_atm(n)
        grid(i,j,k)%y(grid(i,j,k)%num)=yyy_tot_atm(n)
        grid(i,j,k)%z(grid(i,j,k)%num)=zzz_tot_atm(n)
        grid(i,j,k)%ires(grid(i,j,k)%num)=ires_tot_atm(n)
        grid(i,j,k)%anam(grid(i,j,k)%num)=anam_tot_atm(n)
        grid(i,j,k)%rnam(grid(i,j,k)%num)=rnam_tot_atm(n)

      enddo

C ----------------------------------------------------------------------
C now begin a loop seeking to merge regions...

      nrounds=0

666   continue

      nrounds=nrounds+1

C reset counter for regions merged this time

      nmerged=0

C flag that we have looked at no regions this time

      idone(1:num_region)=0

      do n1=1,num_region
        if(idone(n1).eq.1) cycle           ! skip already-done regions
        if(num_mem_region(n1).eq.0) cycle   ! skip now-empty regions
        do n2=n1+1,num_region
          if(idone(n2).eq.1) cycle         ! skip already-done regions
          if(num_mem_region(n2).eq.0) cycle ! skip now-empty regions

C complete-linkage method here - not used but added for future

C measure distances between all members of this pair of regions
C note that we require *all* distances to match the criteria...
C so we skip out as soon as we find any pair that doesn't match

c         imerge=0
c         do m1=1,num_mem_region(n1)
c           do m2=1,num_mem_region(n2)
c             dist2=(xxx_mem_region(n1,m1)-xxx_mem_region(n2,m2))**2+
c    &              (yyy_mem_region(n1,m1)-yyy_mem_region(n2,m2))**2+
c    &              (zzz_mem_region(n1,m1)-zzz_mem_region(n2,m2))**2
c             if(dist2.lt.dist_lo2.or.
c    &           dist2.gt.dist_hi2) then
c               imerge=0
c               goto 777
c             else
c               imerge=1
cc            endif
c           enddo
c         enddo

C single-linkage method here:

          imerge=0
          do m1=1,num_mem_region(n1)
            do m2=1,num_mem_region(n2)
              dist2=(xxx_mem_region(n1,m1)-xxx_mem_region(n2,m2))**2+
     &              (yyy_mem_region(n1,m1)-yyy_mem_region(n2,m2))**2+
     &              (zzz_mem_region(n1,m1)-zzz_mem_region(n2,m2))**2
              if(dist2.ge.dist_lo2.and.
     &           dist2.le.dist_hi2) then
                imerge=1
                goto 777
              endif
            enddo
          enddo

C note that we only get here with imerge=1 if one atom pair matched
C the input distance criteria

777       if(imerge.eq.1) then

C increment #merged this time around

            nmerged=nmerged+1

C copy the members of region n2 in to region n1

            do m2=1,num_mem_region(n2)
            num_mem_region(n1)=num_mem_region(n1)+1
            idr_mem_region(n1,num_mem_region(n1))=idr_mem_region(n2,m2)
            nlg_mem_region(n1,num_mem_region(n1))=nlg_mem_region(n2,m2)
            nam_mem_region(n1,num_mem_region(n1))=nam_mem_region(n2,m2)
            atm_mem_region(n1,num_mem_region(n1))=atm_mem_region(n2,m2)
            xxx_mem_region(n1,num_mem_region(n1))=xxx_mem_region(n2,m2)
            yyy_mem_region(n1,num_mem_region(n1))=yyy_mem_region(n2,m2)
            zzz_mem_region(n1,num_mem_region(n1))=zzz_mem_region(n2,m2)
            bfc_mem_region(n1,num_mem_region(n1))=bfc_mem_region(n2,m2)
            enddo

C now note that we've "done" region n2...

            idone(n2)=1

C and set its #members to zero...

            num_mem_region(n2)=0

          endif
        enddo
      enddo

C now check if we merged any regions this round - if we did, then we
C need to go back and see if any more can be merged

      if(nmerged.gt.0) then
        write(*,*)'in round# ',nrounds,' we merged ',nmerged,' regions'
        goto 666
      endif

C at this point we have merged all possible regions 
C ----------------------------------------------------------------------

C now we may need to add on backbone atoms if required by a ligand - so
C find all backbone atoms, then assign to the nearest region...

      if(iwantbbn.eq.1) then

C first, let's store the original number of members of each region -
C we'll use these to search iteratively - otherwise we'll keep adding
C and adding to the regions...
 
      do n=1,num_region
        num_mem_region_orig(n)=num_mem_region(n)
      enddo

      open(unit=11,file=initial_pdb(1:len1),status="unknown")
90    read(11,'(a80)',end=95)string

C skip non-ATOM entries

      if(string(1:4).ne.'ATOM') goto 90

C now read the entirety of the ATOM entry

      read(string,'(12x,a4,1x,a3,2x,i4,4x,3f8.3,2f6.2)')
     &     anam_tmp,rnam_tmp,rnum_tmp,xt,yt,zt,occ,beta_tmp

C if it's a backbone atom of the appropriate type then assign to nearest
C region if it's within the cutoff...

      if(iwantbbnN.eq.1) then
        if(anam_tmp.eq.' N  ') then
          dist2_best=999999.9
          my_best=-1
          do n=1,num_region
            do m=1,num_mem_region_orig(n)
              dist2=(xt-xxx_mem_region(n,m))**2+
     &              (yt-yyy_mem_region(n,m))**2+
     &              (zt-zzz_mem_region(n,m))**2
              if(dist2.lt.dist2_best) then
                dist2_best=dist2
                my_best=n
              endif
            enddo
          enddo
          if(sqrt(dist2_best).le.bbnadd_dist) then
            write(*,*)'adding BBN N atom to region# ',my_best
            num_mem_region(my_best)=num_mem_region(my_best)+1
            idr_mem_region(my_best,num_mem_region(my_best))=rnum_tmp
            nlg_mem_region(my_best,num_mem_region(my_best))=1
            nam_mem_region(my_best,num_mem_region(my_best))='BBN'
            atm_mem_region(my_best,num_mem_region(my_best))=anam_tmp
            xxx_mem_region(my_best,num_mem_region(my_best))=xt
            yyy_mem_region(my_best,num_mem_region(my_best))=yt
            zzz_mem_region(my_best,num_mem_region(my_best))=zt
            bfc_mem_region(my_best,num_mem_region(my_best))=beta_tmp
          endif
        endif
      elseif(iwantbbnO.eq.1) then
        if(anam_tmp.eq.' O  ') then
          dist2_best=999999.9
          my_best=-1
          do n=1,num_region
            do m=1,num_mem_region_orig(n)
              dist2=(xt-xxx_mem_region(n,m))**2+
     &              (yt-yyy_mem_region(n,m))**2+
     &              (zt-zzz_mem_region(n,m))**2
              if(dist2.lt.dist2_best) then
                dist2_best=dist2
                my_best=n
              endif
            enddo
          enddo
          if(sqrt(dist2_best).le.bbnadd_dist) then
            write(*,*)'adding BBN O atom to region# ',my_best
            num_mem_region(my_best)=num_mem_region(my_best)+1
            idr_mem_region(my_best,num_mem_region(my_best))=rnum_tmp
            nlg_mem_region(my_best,num_mem_region(my_best))=2
            nam_mem_region(my_best,num_mem_region(my_best))='BBN'
            atm_mem_region(my_best,num_mem_region(my_best))=anam_tmp
            xxx_mem_region(my_best,num_mem_region(my_best))=xt
            yyy_mem_region(my_best,num_mem_region(my_best))=yt
            zzz_mem_region(my_best,num_mem_region(my_best))=zt
            bfc_mem_region(my_best,num_mem_region(my_best))=beta_tmp
          endif
        endif
      endif

      goto 90

95    close(11)

      endif

C before going on, we will write out our regions for posterity
C do a simple write to the screen here...

      do n=1,num_region
        do m=1,num_mem_region(n)
          write(*,839)n,m,nam_mem_region(n,m),atm_mem_region(n,m),
     &                    idr_mem_region(n,m),nlg_mem_region(n,m),
     &                    xxx_mem_region(n,m),yyy_mem_region(n,m),
     &                    zzz_mem_region(n,m),bfc_mem_region(n,m)
        enddo
      enddo
839   format('final   region check ',2i6,1x,a3,1x,a4,1x,2i6,4f10.3)

C first, let's figure out if we're likely to have any unregioned CYS
C residues and/or any putative disulfides...

      num_cys_singl=0 ! #CYS that are isolated
      num_cys_doubl=0 ! #CYS in disulfides
      num_cys_lignd=0 ! #CYS with ligands
      num_cys_orfan=0 ! these are evaluated later

      num_disulf=0
      num_region_final=0

      do n1=1,num_region

C skip now-empty regions

        if(num_mem_region(n1).eq.0) cycle

C go on and increment the # of regions

        num_region_final=num_region_final+1

C look for single cysteines

        if(num_mem_region(n1).eq.1.and.
     &     nam_mem_region(n1,1).eq.'CYS') 
     &     num_cys_singl=num_cys_singl+1

C evaluate potential disulfides - note that some of these may well
C eventually end up as liganded to a metal cluster

        do m1=1,num_mem_region(n1)
          if(nam_mem_region(n1,m1).ne.'CYS') cycle
          do m2=m1+1,num_mem_region(n1)
            if(nam_mem_region(n1,m2).ne.'CYS') cycle
            dist2=(xxx_mem_region(n1,m1)-xxx_mem_region(n1,m2))**2+
     &            (yyy_mem_region(n1,m1)-yyy_mem_region(n1,m2))**2+ 
     &            (zzz_mem_region(n1,m1)-zzz_mem_region(n1,m2))**2
            if(dist2.le.disulf_dist2) then
              num_disulf=num_disulf+1
              ncl_disulf(num_disulf)=num_region_final
              idr_disulf(num_disulf)=m1
              jdr_disulf(num_disulf)=m2
              dis_disulf(num_disulf)=sqrt(dist2)
              num_cys_doubl=num_cys_doubl+2
            endif
          enddo
        enddo
      enddo

C write them out to file outfile6

      open(unit=32,file=outfile6,status="unknown")

      do n=1,num_disulf
        write(32,551)initial_pdb(1:30),n,num_disulf,
     &               ncl_disulf(n),
     &               idr_disulf(n),
     &               jdr_disulf(n),
     &               dis_disulf(n)
551   format('pdbname ',a30,' disulf# ',i6,' out_of ',i6,
     &      ' disulfides; is_in_region# ',i6,
     &      ' and_involves_res_nums: ',2i6,' distance: ',f10.3)
      enddo

      close(32)

C now also write out all region information to file outfile5

      open(unit=32,file=outfile5,status="unknown")
      open(unit=33,file=outfile7,status="unknown")

      num_region_final=0
      do n1=1,num_region

C skip now-empty regions

        if(num_mem_region(n1).eq.0) cycle

C go on and increment the # of regions

        num_region_final=num_region_final+1
        idr_last=-1 
        num_mem_final=0

C here we do a loop to find out the identities of all unique *residues*
C in the region (remember that HIS residues might have two atoms from
C the same residue in the same region)

        do m1=1,num_mem_region(n1)
          if(idr_mem_region(n1,m1).ne.idr_last) then
            num_mem_final=num_mem_final+1
            idr_mem_final(num_mem_final)=idr_mem_region(n1,m1)
            nam_mem_final(num_mem_final)=nam_mem_region(n1,m1)
            idr_last=idr_mem_region(n1,m1)
          endif
        enddo

C now write out to the screen - note that we write out the total number
C of members (atoms) in the region (num_mem_region(n1)) *and* the total
C number of unique residues in the region (num_mem_final)

        write(*,101)initial_pdb(1:30),num_region_final,
     &              num_mem_region(n1),
     &              num_mem_final,
     &             (nam_mem_final(j),idr_mem_final(j),j=1,num_mem_final)

        write(32,101)initial_pdb(1:30),num_region_final,
     &              num_mem_region(n1), ! # of atoms in region
     &              num_mem_final,     ! # of residues in region
     &             (nam_mem_final(j),idr_mem_final(j),j=1,num_mem_final)

C if this region contains only 4 CYS and 4 CYS only, then we will check
C to see if all distances are within the cutoff and, if so, record it

        ispecial_region(n1)=0

        if(num_mem_final.eq.4) then

C default to assuming that this *is* a special region

          ispecial_region(n1)=1

C first skip it if any are not CYS

          do m1=1,num_mem_region(n1)
            if(nam_mem_region(n1,m1).ne.'CYS') then
              ispecial_region(n1)=0
              goto 727
            endif
          enddo

C measure distances between all members of this region - if any of them
C are outside the range then set ispecial to zero and exit 

          do m1=1,num_mem_region(n1)
            do m2=m1+1,num_mem_region(n1)
              dist2=(xxx_mem_region(n1,m1)-xxx_mem_region(n1,m2))**2+
     &              (yyy_mem_region(n1,m1)-yyy_mem_region(n1,m2))**2+
     &              (zzz_mem_region(n1,m1)-zzz_mem_region(n1,m2))**2
              if(dist2.lt.dist_lo2.or.
     &           dist2.gt.dist_hi2) then
                ispecial_region(n1)=0
                goto 727
              endif
            enddo
          enddo

727       continue

        endif 
          
      enddo

101   format('pdbname ',a30,' FINAL REGION INFO ;',
     &      ' region# ',i6,' contains: ',i6,' atoms & ',i6,
     &      ' residues: ',
     &        100(a3,i6,' | '))

      close(32)

C ----------------------------------------------------------------------
C now we ask which ligand type fits best into each region
C ----------------------------------------------------------------------

      num_region_final=0
      do n1=1,num_region

C skip now-empty regions

        if(num_mem_region(n1).eq.0) cycle

C go on and increment the # of regions

        num_region_final=num_region_final+1

C skip regions that have two few members to match any ligand
C note that before skipping we ask if there are any orfan CYS residues
C in those regions that contain >1 member (CYS residues in regions
C with 1 member are added to num_cys_singl instead)

        if(num_mem_region(n1).lt.nmtc_min) then
          if(num_mem_region(n1).eq.1) cycle
          do m1=1,num_mem_region(n1)
            if(nam_mem_region(n1,m1).eq.'CYS') then
              num_cys_orfan=num_cys_orfan+1
            endif
          enddo
          cycle
        endif

C we're not done yet! we also need to skip regions that, while they may
C have enough members, don't have enough unique residue numbers to match

        num_mem_final=0
        do m1=1,num_mem_region(n1)
          iaccept=1
          do k1=1,num_mem_final
            if(idr_mem_region(n1,m1).eq.idr_mem_final(k1)) then
              iaccept=0
              exit
            endif
          enddo
          if(iaccept.eq.1) then
            num_mem_final=num_mem_final+1
            idr_mem_final(num_mem_final)=idr_mem_region(n1,m1)
          endif
        enddo

C again, before quitting on this region we increment num_cys_orfan...

        if(num_mem_final.lt.nmtc_min) then
          if(num_mem_region(n1).eq.1) cycle
          do m1=1,num_mem_region(n1)
            if(nam_mem_region(n1,m1).eq.'CYS') then
              num_cys_orfan=num_cys_orfan+1
            endif
          enddo
          cycle
          cycle
        endif

C otherwise go on

        rms_fin=999999.9
        write(*,*)'WORKING on region ',num_region_final,
     &           num_mem_region(n1),
     &          (nam_mem_region(n1,m1),m1=1,num_mem_region(n1))

C first, record all members as unassigned (i.e. not "dun")
C note that dun_mem_region is a temporary array reset for each region
C same with mlg_mem_region which records the # of different ligands
C attached to this atom and ilg_mem_region which records the ligsite#

        dun_mem_region(1:num_mem_region(n1))=0
        mlg_mem_region(1:num_mem_region(n1))=0
        ilg_mem_region(1:num_mem_region(n1),1:3)=0

C also set # of iterations to zero...

        num_iters=0

C come back here if some members of the region still unassigned

444     continue

        num_iters=num_iters+1

C loop over all ligand types

        do nnn=1,nlig

C skip this ligand type if it has too many match atoms - it can't
C possibly be a good match for this particular region

          if(nmtc(nnn).gt.num_mem_region(n1)) then
            write(*,*)'skipping ligand type ',nnn
            rms_best(nnn)=999999.9
            cycle
          endif

C otherwise set default values for rms_best,i_best etc

          ifitbutclashed=0
          rms_best(nnn)=999999.9
          i_best(nnn)=-1
          j_best(nnn)=-1
          k_best(nnn)=-1
          l_best(nnn)=-1

C get the coords for the (moving) ligand - all atoms

          do mmm=1,nmov(nnn)
            xm(mmm)=xlig(nnn,mmm)
            ym(mmm)=ylig(nnn,mmm)
            zm(mmm)=zlig(nnn,mmm)
          enddo

C now loop over all possible combinations of fixed atoms in this region

          num_comb=0

          do i=1,num_mem_region(n1)

C skip if this combo has been "dun"

ccc         if(dun_mem_region(i).ge.1) cycle
            if(mlg_mem_region(i).ge.nlg_mem_region(n1,i)) cycle

            do j=1,num_mem_region(n1)      

C note that we skip combos where the same atom has been selected twice

              if(j.eq.i) cycle

C *and* we skip ones where two atoms are from the same residue

              if(idr_mem_region(n1,j).eq.idr_mem_region(n1,i)) cycle

C *and* we skip if this combo has been "dun"

ccc           if(dun_mem_region(j).ge.1) cycle
              if(mlg_mem_region(j).ge.nlg_mem_region(n1,j)) cycle

              do k=1,num_mem_region(n1)      

                if(k.eq.i) cycle
                if(k.eq.j) cycle
                if(idr_mem_region(n1,k).eq.idr_mem_region(n1,i)) cycle
                if(idr_mem_region(n1,k).eq.idr_mem_region(n1,j)) cycle
ccc             if(dun_mem_region(k).ge.1) cycle
                if(mlg_mem_region(k).ge.nlg_mem_region(n1,k)) cycle


C we need to have two possible paths here depending on whether
C #match atoms (num_mtc) is 3 or 4 - the default case is 4, so we need
C i,j,k,l - but for 3 we only need i,j,k (we need to try at least 1
C value of l though to make sure we go through the loop)

                if(nmtc(nnn).eq.3) then
                  lbeg=1
                  lend=1
                elseif(nmtc(nnn).ge.4) then
                  lbeg=1
                  lend=num_mem_region(n1)
                endif

                do l=lbeg,lend

C we only check for repeated entries when num_mem_region>=4 (see above)
C same thing applies if we check that this position has been "dun"

                if(nmtc(nnn).ge.4) then
                  if(l.eq.i) cycle
                  if(l.eq.j) cycle
                  if(l.eq.k) cycle
                  if(idr_mem_region(n1,l).eq.idr_mem_region(n1,i)) cycle
                  if(idr_mem_region(n1,l).eq.idr_mem_region(n1,j)) cycle
                  if(idr_mem_region(n1,l).eq.idr_mem_region(n1,k)) cycle
ccc               if(dun_mem_region(l).ge.1) cycle
                  if(mlg_mem_region(l).ge.nlg_mem_region(n1,l)) cycle
                endif

C we need to have two possible paths here depending on whether
C #match atoms (num_mtc) is 3 or 4 - the default case is 4, so we need
C i,j,k,l - but for 3 we only need i,j,k (we need to try at least 1
C value of l though to make sure we go through the loop)

                if(nmtc(nnn).le.4) then
                  mbeg=1
                  mend=1
                elseif(nmtc(nnn).ge.5) then
                  mbeg=1
                  mend=num_mem_region(n1)
                endif

                do m=mbeg,mend

C we only check for repeated entries when num_mem_region>=4 (see above)
C same thing applies if we check that this position has been "dun"

                if(nmtc(nnn).ge.5) then
                  if(m.eq.i) cycle
                  if(m.eq.j) cycle
                  if(m.eq.k) cycle
                  if(m.eq.l) cycle
                  if(idr_mem_region(n1,m).eq.idr_mem_region(n1,i)) cycle
                  if(idr_mem_region(n1,m).eq.idr_mem_region(n1,j)) cycle
                  if(idr_mem_region(n1,m).eq.idr_mem_region(n1,k)) cycle
                  if(idr_mem_region(n1,m).eq.idr_mem_region(n1,l)) cycle
ccc               if(dun_mem_region(m).ge.1) cycle
                  if(mlg_mem_region(m).ge.nlg_mem_region(n1,m)) cycle
                endif

C we need to have two possible paths here depending on whether
C #match atoms (num_mtc) is 3 or 4 - the default case is 4, so we need
C i,j,k,l - but for 3 we only need i,j,k (we need to try at least 1
C value of l though to make sure we go through the loop)

                if(nmtc(nnn).le.5) then
                  nbeg=1
                  nend=1
                elseif(nmtc(nnn).ge.6) then
                  nbeg=1
                  nend=num_mem_region(n1)
                endif

                do n=nbeg,nend

C we only check for repeated entries when num_mem_region>=4 (see above)
C same thing applies if we check that this position has been "dun"

                if(nmtc(nnn).ge.6) then
                  if(n.eq.i) cycle
                  if(n.eq.j) cycle
                  if(n.eq.k) cycle
                  if(n.eq.l) cycle
                  if(n.eq.m) cycle
                  if(idr_mem_region(n1,n).eq.idr_mem_region(n1,i)) cycle
                  if(idr_mem_region(n1,n).eq.idr_mem_region(n1,j)) cycle
                  if(idr_mem_region(n1,n).eq.idr_mem_region(n1,k)) cycle
                  if(idr_mem_region(n1,n).eq.idr_mem_region(n1,l)) cycle
                  if(idr_mem_region(n1,n).eq.idr_mem_region(n1,m)) cycle
ccc               if(dun_mem_region(n).ge.1) cycle
                  if(mlg_mem_region(n).ge.nlg_mem_region(n1,n)) cycle
                endif

C finally, we also need to check that we have the right number of
C CYS and HIS atoms for this particular combination... - note that this
C code relies on us already having filtered out those cases where two
C atoms of the same residue are present in the current combo

                  ncys_tmp=0
                  nhis_tmp=0
                  nasp_tmp=0
                  ndoe_tmp=0
                  nbbn_tmp=0
                  if(nam_mem_region(n1,i).eq.'CYS') ncys_tmp=ncys_tmp+1
                  if(nam_mem_region(n1,j).eq.'CYS') ncys_tmp=ncys_tmp+1
                  if(nam_mem_region(n1,k).eq.'CYS') ncys_tmp=ncys_tmp+1
                  if(nam_mem_region(n1,i).eq.'HIS') nhis_tmp=nhis_tmp+1
                  if(nam_mem_region(n1,j).eq.'HIS') nhis_tmp=nhis_tmp+1
                  if(nam_mem_region(n1,k).eq.'HIS') nhis_tmp=nhis_tmp+1
                  if(nam_mem_region(n1,i).eq.'ASP') nasp_tmp=nasp_tmp+1
                  if(nam_mem_region(n1,j).eq.'ASP') nasp_tmp=nasp_tmp+1
                  if(nam_mem_region(n1,k).eq.'ASP') nasp_tmp=nasp_tmp+1
                  if(nam_mem_region(n1,i).eq.'DOE') ndoe_tmp=ndoe_tmp+1
                  if(nam_mem_region(n1,j).eq.'DOE') ndoe_tmp=ndoe_tmp+1
                  if(nam_mem_region(n1,k).eq.'DOE') ndoe_tmp=ndoe_tmp+1
                  if(nam_mem_region(n1,i).eq.'BBN') nbbn_tmp=nbbn_tmp+1
                  if(nam_mem_region(n1,j).eq.'BBN') nbbn_tmp=nbbn_tmp+1
                  if(nam_mem_region(n1,k).eq.'BBN') nbbn_tmp=nbbn_tmp+1
                  if(nmtc(nnn).ge.4) then
                  if(nam_mem_region(n1,l).eq.'CYS') ncys_tmp=ncys_tmp+1
                  if(nam_mem_region(n1,l).eq.'HIS') nhis_tmp=nhis_tmp+1
                  if(nam_mem_region(n1,l).eq.'ASP') nasp_tmp=nasp_tmp+1
                  if(nam_mem_region(n1,l).eq.'DOE') ndoe_tmp=ndoe_tmp+1
                  if(nam_mem_region(n1,l).eq.'BBN') nbbn_tmp=nbbn_tmp+1
                  endif
                  if(nmtc(nnn).ge.5) then
                  if(nam_mem_region(n1,m).eq.'CYS') ncys_tmp=ncys_tmp+1
                  if(nam_mem_region(n1,m).eq.'HIS') nhis_tmp=nhis_tmp+1
                  if(nam_mem_region(n1,m).eq.'ASP') nasp_tmp=nasp_tmp+1
                  if(nam_mem_region(n1,m).eq.'DOE') ndoe_tmp=ndoe_tmp+1
                  if(nam_mem_region(n1,m).eq.'BBN') nbbn_tmp=nbbn_tmp+1
                  endif
                  if(nmtc(nnn).ge.6) then
                  if(nam_mem_region(n1,n).eq.'CYS') ncys_tmp=ncys_tmp+1
                  if(nam_mem_region(n1,n).eq.'HIS') nhis_tmp=nhis_tmp+1
                  if(nam_mem_region(n1,n).eq.'ASP') nasp_tmp=nasp_tmp+1
                  if(nam_mem_region(n1,n).eq.'DOE') ndoe_tmp=ndoe_tmp+1
                  if(nam_mem_region(n1,n).eq.'BBN') nbbn_tmp=nbbn_tmp+1
                  endif

                  if(ncys_tmp.ne.ncys_lig(nnn).or.
     &               nhis_tmp.ne.nhis_lig(nnn).or.
     &               nasp_tmp.ne.nasp_lig(nnn).or.
     &               ndoe_tmp.ne.ndoe_lig(nnn).or.
     &               nbbn_tmp.ne.nbbn_lig(nnn)) then
c                   write(*,*)'skipping combo LIG ',
c    &                        ncys_lig(nnn),nhis_lig(nnn),' TMP ',
c    &                        ncys_tmp,nhis_tmp
                    cycle
                  endif

C otherwise, let's go on and do a superposition
C first, set the coords of the fixed atoms
       
                  num_comb=num_comb+1
                  xf(1)=xxx_mem_region(n1,i)
                  yf(1)=yyy_mem_region(n1,i)
                  zf(1)=zzz_mem_region(n1,i)
                  xf(2)=xxx_mem_region(n1,j)
                  yf(2)=yyy_mem_region(n1,j)
                  zf(2)=zzz_mem_region(n1,j)
                  xf(3)=xxx_mem_region(n1,k)
                  yf(3)=yyy_mem_region(n1,k)
                  zf(3)=zzz_mem_region(n1,k)
                  if(nmtc(nnn).ge.4) then
                    xf(4)=xxx_mem_region(n1,l)
                    yf(4)=yyy_mem_region(n1,l)
                    zf(4)=zzz_mem_region(n1,l)
                  endif
                  if(nmtc(nnn).ge.5) then
                    xf(5)=xxx_mem_region(n1,m)
                    yf(5)=yyy_mem_region(n1,m)
                    zf(5)=zzz_mem_region(n1,m)
                  endif
                  if(nmtc(nnn).ge.6) then
                    xf(6)=xxx_mem_region(n1,n)
                    yf(6)=yyy_mem_region(n1,n)
                    zf(6)=zzz_mem_region(n1,n)
                  endif
               
C now superimpose

                  call pdbsup(nmov(nnn),nmtc(nnn),xm,ym,zm,
     &                        xf,yf,zf,rms)
c                 write(*,*)'CHECK RMS ',nnn,rms

C now do a check to make sure we have no steric clashes
C note that we only do this if the rms<=rms_thresh...

                  if(rms.le.rms_thresh(nnn)) then

                    iclash=0

C note that we loop over the non-match atoms only...

                    dist2_min=999999.9

                    do mmm=nmtc(nnn)+1,nmov(nnn)

                      ixx=int((xm(mmm)-xmin)/rcll)+1
                      jxx=int((ym(mmm)-ymin)/rcll)+1
                      kxx=int((zm(mmm)-zmin)/rcll)+1
                      ilo=ixx-1
                      jlo=jxx-1
                      klo=kxx-1
                      ihi=ixx+1
                      jhi=jxx+1
                      khi=kxx+1
 
                      if(ilo.lt.1) ilo=1
                      if(jlo.lt.1) jlo=1
                      if(klo.lt.1) klo=1
                      if(ihi.gt.1) ihi=numx
                      if(jhi.gt.1) jhi=numy
                      if(khi.gt.1) khi=numz

                      do ixx=ilo,ihi
                      do jxx=jlo,jhi
                      do kxx=klo,khi
                      do nop=1,grid(ixx,jxx,kxx)%num

C check whether the current grid atom is one of the atoms used to place
C the ligand - if it is, then we should skip it...

                        rnum_tmp=grid(ixx,jxx,kxx)%ires(nop)
                        anam_tmp=grid(ixx,jxx,kxx)%anam(nop)
                        rnam_tmp=grid(ixx,jxx,kxx)%rnam(nop)

                        if(rnum_tmp.eq.idr_mem_region(n1,i).and.
     &                     anam_tmp.eq.atm_mem_region(n1,i).and.
     &                     rnam_tmp.eq.nam_mem_region(n1,i)) cycle
                        if(rnum_tmp.eq.idr_mem_region(n1,j).and.
     &                     anam_tmp.eq.atm_mem_region(n1,j).and.
     &                     rnam_tmp.eq.nam_mem_region(n1,j)) cycle
                        if(rnum_tmp.eq.idr_mem_region(n1,k).and.
     &                     anam_tmp.eq.atm_mem_region(n1,k).and.
     &                     rnam_tmp.eq.nam_mem_region(n1,k)) cycle
                        if(nmtc(nnn).ge.4) then
                          if(rnum_tmp.eq.idr_mem_region(n1,l).and.
     &                       anam_tmp.eq.atm_mem_region(n1,l).and.
     &                       rnam_tmp.eq.nam_mem_region(n1,l)) cycle
                        endif
                        if(nmtc(nnn).ge.5) then
                          if(rnum_tmp.eq.idr_mem_region(n1,m).and.
     &                       anam_tmp.eq.atm_mem_region(n1,m).and.
     &                       rnam_tmp.eq.nam_mem_region(n1,m)) cycle
                        endif
                        if(nmtc(nnn).ge.6) then
                          if(rnum_tmp.eq.idr_mem_region(n1,n).and.
     &                       anam_tmp.eq.atm_mem_region(n1,n).and.
     &                       rnam_tmp.eq.nam_mem_region(n1,n)) cycle
                        endif

C if we get here then we need to measure the distance for clash check

                        dist2=(xm(mmm)-grid(ixx,jxx,kxx)%x(nop))**2+
     &                        (ym(mmm)-grid(ixx,jxx,kxx)%y(nop))**2+
     &                        (zm(mmm)-grid(ixx,jxx,kxx)%z(nop))**2
                        dist2_min=min(dist2_min,dist2)
                        if(dist2.lt.clash2) then
                          iclash=1
                          goto 222
                        endif

                      enddo
                      enddo
                      enddo
                      enddo

                    enddo

C now check against previously-added ligand atoms - note that we only
C get here if there was no clash with protein - note also that we mark
C these clashes as iclash=2...

                    dist2_min=999999.9

                    do mmm=nmtc(nnn)+1,nmov(nnn)
                      do mm2=1,nfin_atm
                        dist2=(xfin(mm2)-xm(mmm))**2+
     &                        (yfin(mm2)-ym(mmm))**2+
     &                        (zfin(mm2)-zm(mmm))**2
                        dist2_min=min(dist2_min,dist2)
                        if(dist2.lt.clash_lig_dist2) then
                          iclash=2
                          goto 222
                        endif
                      enddo
                    enddo

C come here if we jumped out with any kind of clash...

222                 continue
  
                    if(iclash.eq.1) then

C note that we only set ifitbutclashed to 1 if it has NOT already had a
C good match that didn't clash... - otherwise we stand a chance of
C thinking that a *later* rmsd-match is a clasher and that this
C overrides the fact that an *earlier* rmsd-match fit great

                      if(ifitbutclashed.ne.-1) then
                        ifitbutclashed=1
                      endif
ccc                   write(*,*)'Clash Detected! ',rms,sqrt(dist2)
                      cycle
                    elseif(iclash.eq.2) then
ccc                   write(*,*)'Clash Detected! ',rms,sqrt(dist2)
                      cycle
                    endif

                  endif

C if the superposition is the best yet for this ligand type *and* it has
C not steric clashes then update all the information for it...
C remember that we only get here if iclash=0...

                  if(rms.lt.rms_best(nnn)) then

C note that we record here that we found at least one case where the
C rmsd was good and the steric fit was fine too

                    ifitbutclashed=-1
                    rms_best(nnn)=rms
                    i_best(nnn)=i
                    j_best(nnn)=j
                    k_best(nnn)=k
                    l_best(nnn)=l
                    m_best(nnn)=m
                    n_best(nnn)=n
                  endif

                enddo
                enddo
                enddo
              enddo
            enddo
          enddo

C write out the best value to the screen for checking - but only do this
C for ligands that were actually tested...
C note that we also write out whether the ligand sterically clashes with
C the rest of the protein - we do NOT count if the clash is with a
C previously-added ligand...

          if(num_comb.ge.1) then
            if(ifitbutclashed.eq.-1) then
              write(*,611)initial_pdb(1:30),
     &                    num_region_final,nnn,rms_best(nnn)
611           format('pdbname ',a30,' region# ',i6,' ligand# ',i6,
     &                ' rms_best = ',f12.3,' steric_fits_great')
            elseif(ifitbutclashed.eq.1) then
              write(*,711)initial_pdb(1:30),
     &                    num_region_final,nnn,rms_best(nnn)
711           format('pdbname ',a30,' region# ',i6,' ligand# ',i6,
     &                ' rms_best = ',f12.3,' but_clashes_with_protein')
            endif
          endif

        enddo

C NEW APPROACH

C now that all ligands have had their best shot we want to pick a final
C answer - we'll filter out those ligands whose best rmsd is above the
C threshold and then we'll keep the FIRST type of ligand that matches
C the criterion - easy to implement but needs testing

        nmtc_best=0
        nnn_fin=0

        do nnn=1,nlig

          if(rms_best(nnn).gt.rms_thresh(nnn)) then
            cycle
          else
            write(*,613)n1,nnn,rms_best(nnn),rms_thresh(nnn)
613         format('for region# ',i6,' new best ligand is ',i6,
     &            ' with rmsd of ',f8.3,' beats thresh ',f8.3)
            nmtc_best=nmov(nnn)
            nnn_fin=nnn
            rms_fin=rms_best(nnn)
            exit
          endif

        enddo

C need to quit completely if no good match - if there is a match we will
C have set nmtc_best in the lines above...

        if(nmtc_best.eq.0) then
          num_orf=0
          do m=1,num_mem_region(n1)
            if(dun_mem_region(m).eq.0) then
              num_orf=num_orf+1
              if(nam_mem_region(n1,m).eq.'CYS') then
                num_cys_orfan=num_cys_orfan+1
              endif
            endif
          enddo
          write(*,*)'no more matches for region# ',num_region_final,
     &             ' num_orfs = ',num_orf
          goto 333
        endif

C now we can now write out the transformed coordinates of the best-fit
C ligand to a pdb file and we can determine other features of interest

C first need to restore the i,j,k,l of the best fit

        if(rms_fin.le.rms_thresh(nnn_fin)) then

          write(*,615)initial_pdb(1:30),num_region_final,
     &                num_iters,nnn_fin,rms_fin
615       format('pdbname ',a30,
     &        ' for_region# ',i6,' and_iteration_# ',i6,
     &        ' FINAL_ANSWER_is ',i6,
     &        ' with_RMSD_of ',f8.3)

          nligsites=nligsites+1
          if(nligsites.eq.1) then
            open(unit=13,file=outfile3,status='unknown')
          endif

          nnn=nnn_fin
          i=i_best(nnn)
          j=j_best(nnn)
          k=k_best(nnn)
          if(nmtc(nnn).ge.4) then
            l=l_best(nnn)
          endif
          if(nmtc(nnn).ge.5) then
            m=m_best(nnn)
          endif
          if(nmtc(nnn).ge.6) then
            n=n_best(nnn)
          endif

C let's flag that these positions have been taken by a region
C note that we assign them the number of ligands found thus far

          mlg_mem_region(i)=mlg_mem_region(i)+1
          mlg_mem_region(j)=mlg_mem_region(j)+1
          mlg_mem_region(k)=mlg_mem_region(k)+1
          ilg_mem_region(i,mlg_mem_region(i))=nligsites
          ilg_mem_region(j,mlg_mem_region(j))=nligsites
          ilg_mem_region(k,mlg_mem_region(k))=nligsites
       
          dun_mem_region(i)=nligsites
          dun_mem_region(j)=nligsites
          dun_mem_region(k)=nligsites

C crude flag to see if code is working...

          if(mlg_mem_region(i).eq.nlg_mem_region(n1,i).or.
     &       mlg_mem_region(j).eq.nlg_mem_region(n1,j).or.
     &       mlg_mem_region(k).eq.nlg_mem_region(n1,k)) then
            write(*,943)initial_pdb(1:30)
943         format('FOUND SHARED SITE IN ',a30)
          endif

          if(nmtc(nnn).ge.4) then
            dun_mem_region(l)=nligsites
            mlg_mem_region(l)=mlg_mem_region(l)+1
            ilg_mem_region(l,mlg_mem_region(l))=nligsites
          endif
          if(nmtc(nnn).ge.5) then
            dun_mem_region(m)=nligsites
            mlg_mem_region(m)=mlg_mem_region(m)+1
            ilg_mem_region(m,mlg_mem_region(m))=nligsites
          endif
          if(nmtc(nnn).ge.6) then
            dun_mem_region(n)=nligsites
            mlg_mem_region(n)=mlg_mem_region(n)+1
            ilg_mem_region(n,mlg_mem_region(n))=nligsites
          endif

C and let's figure out if we're done (use this later) - note that here
C we can still use dun_mem_region as this will record whether at least
C one ligand has been attached to this site...

          iam_not_dun=0
          do mx=1,num_mem_region(n1)
            if(dun_mem_region(mx).eq.0) then
              iam_not_dun=1
            else
              if(nam_mem_region(n1,mx).eq.'CYS') then
                num_cys_lignd=num_cys_lignd+1
              endif
            endif
          enddo 

C let's increment ifound and nfound for this ligand type

          if(ifound(nnn).eq.0) ifound(nnn)=1
          nfound(nnn)=nfound(nnn)+1

C let's get the min, max, average B-factor of the residues used in the
C fit - we might find later a relationship between this and the rmsd...

          bfc_max=-999.9
          bfc_min= 999.9
          bfc_ave=   0.0
          if(nmtc(nnn).ge.3) then
            bfc_max=max(bfc_max,bfc_mem_region(n1,i))
            bfc_max=max(bfc_max,bfc_mem_region(n1,j))
            bfc_max=max(bfc_max,bfc_mem_region(n1,k))
            bfc_min=min(bfc_min,bfc_mem_region(n1,i))
            bfc_min=min(bfc_min,bfc_mem_region(n1,j))
            bfc_min=min(bfc_min,bfc_mem_region(n1,k))
            bfc_ave=bfc_ave+bfc_mem_region(n1,i)
            bfc_ave=bfc_ave+bfc_mem_region(n1,j)
            bfc_ave=bfc_ave+bfc_mem_region(n1,k)
          endif
          if(nmtc(nnn).ge.4) then
            bfc_max=max(bfc_max,bfc_mem_region(n1,l))
            bfc_min=min(bfc_min,bfc_mem_region(n1,l))
            bfc_ave=bfc_ave+bfc_mem_region(n1,l)
          endif
          if(nmtc(nnn).ge.5) then
            bfc_max=max(bfc_max,bfc_mem_region(n1,m))
            bfc_min=min(bfc_min,bfc_mem_region(n1,m))
            bfc_ave=bfc_ave+bfc_mem_region(n1,m)
          endif
          if(nmtc(nnn).ge.6) then
            bfc_max=max(bfc_max,bfc_mem_region(n1,n))
            bfc_min=min(bfc_min,bfc_mem_region(n1,n))
            bfc_ave=bfc_ave+bfc_mem_region(n1,n)
          endif
          bfc_ave=bfc_ave/real(nmtc(nnn))

C for writing out put the res#s into char6 format
C note that we also take care to only pull out the sites that are part
C of the current ligand's binding site...

          num_mem_tmp=0
          do mx=1,num_mem_region(n1)

C if now ligands attached to this site then skip it

            if(mlg_mem_region(mx).eq.0) cycle

C otherwise we need to make sure that we're using the most recent ligand
C attached to this particular site...

ccc         if(dun_mem_region(mx).ne.nligsites) cycle
            if(mlg_mem_region(mx).gt.0) then
              if(ilg_mem_region(mx,mlg_mem_region(mx)).ne.
     &           nligsites) cycle
            endif

            num_mem_tmp=num_mem_tmp+1
            nam_mem_tmp(num_mem_tmp)=nam_mem_region(n1,mx)
            write(char6,'(i6)')idr_mem_region(n1,mx)
            do m2=1,6
              if(char6(m2:m2).eq.' ') char6(m2:m2)='0'
            enddo
            char6_arr(num_mem_tmp)=char6
          enddo

C let's also get the rmsd of the best fit and the second best fit...
C note that we only compare our best fit with other ligands that have
C the same number of match atoms: it doesn't make sense to compare a
C 4-atom fit with a 3-atom fit as the latter will always be better

          delta_rms=99.999
          do nn2=1,nlig

            if(nn2.eq.nnn_fin) cycle        ! skip the best one
            if(nmtc(nn2).lt.nmtc(nnn_fin)) cycle ! skip ones with   
                                                 ! fewer match atoms
            if(rms_best(nn2).gt.99.9) cycle ! skip ones with no rms

            delta=rms_best(nn2)-rms_fin
            if(delta.lt.delta_rms) delta_rms=delta

          enddo

          if(delta_rms.lt.99.999) then
            rms_2fin=rms_fin+delta_rms
          else
            rms_2fin=99.999
          endif

C finally, if this is a special region (i.e. 4CYS) let's get all RMSDs
C and write them out to file outfile7
C note that we need to be very careful in cases where ths special
C region did NOT match to one of our 4CYS ligands - in such cases it's
C possible that it'll have matched to a 3CYS ligand - so we'll have
C separate code to handle such cases

          if(ispecial_region(n1).eq.1) then
            write(*,*)
            write(*,*)'SPECIAL REGION ',num_region_final
            write(*,*)
            ntyps=0

C if a 4CYS ligand did match the rms_thresh then everything's easy:

            if(nmtc(nnn_fin).eq.4) then

              do nn2=1,nlig
                if(nmtc(nn2).lt.nmtc(nnn_fin)) cycle ! skip ones with   
                                                     ! fewer match atoms
                if(rms_best(nn2).gt.99.9) cycle ! skip ones with no rms
                ntyps=ntyps+1
                rms_typ(ntyps)=rms_best(nn2)
              enddo

C note that we assume for formating that there will be only THREE
C ligands and there will be FOUR CYSs - TEMP TEMP TEMP

              write(33,352)initial_pdb(1:30),num_region_final,
     &             (nam_mem_tmp(mx),char6_arr(mx),mx=1,4),
     &             (rms_typ(lx),lx=1,ntyps),
     &             (bfc_mem_region(n1,l1),l1=1,num_mem_region(n1))

352           format('pdbname ',a30,
     &        ' 4CYSonly region# ',i6,
     &        ' SUCCESS ',
     &        ' has_members: ',4(a3,'_',a6,' | '),
     &        ' and_best_ligand_RMSDs: ',3f7.3,
     &        ' and_bfactors: ',4f7.3)

C if a 4CYS ligand did NOT match the rms_thresh then really all we want
C to do here is write out their RMSDs (and get their resnumbers correct)

            elseif(nmtc(nnn_fin).lt.4) then

              do nn2=1,nlig
                if(nmtc(nn2).lt.4) cycle        ! skip ones with   
                                                ! fewer match atoms
                if(rms_best(nn2).gt.99.9) cycle ! skip ones with no rms
                ntyps=ntyps+1
                rms_typ(ntyps)=rms_best(nn2)
              enddo

C now we also need to make sure that the res#s are written correctly

              do mx=1,num_mem_region(n1)
                write(char6,'(i6)')idr_mem_region(n1,mx)
                do m2=1,6
                  if(char6(m2:m2).eq.' ') char6(m2:m2)='0'
                enddo
                char6_brr(mx)=char6
              enddo

C note that we assume for formating that there will be only THREE
C ligands and there will be FOUR CYSs - TEMP TEMP TEMP

              write(33,353)initial_pdb(1:30),num_region_final,
     &             (nam_mem_region(n1,mx),char6_brr(mx),mx=1,4),
     &             (rms_typ(lx),lx=1,ntyps),
     &             (bfc_mem_region(n1,l1),l1=1,num_mem_region(n1))

353           format('pdbname ',a30,
     &        ' 4CYSonly region# ',i6,
     &        ' FAILURE ',
     &        ' has_members: ',4(a3,'_',a6,' | '),
     &        ' and_best_ligand_RMSDs: ',3f7.3,
     &        ' and_bfactors: ',4f7.3)

            endif

          endif

C write new ligand and region to info file

          write(31,351)initial_pdb(1:30),nligsites,nnn,
     &            ligand_pdb(nnn)(1:19),n1,num_iters,
     &            bfc_min,bfc_ave,bfc_max,
     &            rms_fin,rms_2fin,delta_rms,
     &           (nam_mem_tmp(mx),char6_arr(mx),
     &        mx=1,num_mem_tmp)

351   format('pdbname ',a30,
     &      ' FINAL ligandsite# ',i6,' IS TYPE:  ',i6,
     &      ' ligand_name ',a19,
     &      ' region# ',i6,' iteration# ',i6,
     &      ' bfc_min,bfc_ave,bfc_max: ',3f6.2,
     &      ' rms_best,rms_next_best,delta_rms: ',3f7.3,' | ',
     &        100(a3,'_',a6,' | '))

          write(*,301)initial_pdb(1:30),num_region_final,nnn,nligsites,
     &            bfc_min,bfc_ave,bfc_max,
     &            rms_fin,rms_2fin,delta_rms,
     &           (nam_mem_tmp(mx),char6_arr(mx),
ccc  &           (nam_mem_region(n1,mx),char6_arr(mx),
     &        mx=1,num_mem_tmp)
ccc  &        mx=1,num_mem_region(n1))

301   format('pdbname ',a30,
     &      ' FINAL    region# ',i6,' IS TYPE:  ',i6,
     &      ' and_is_ligsite# ',i6,
     &      ' bfc_min,bfc_ave,bfc_max: ',3f6.2,
     &      ' rms_best,rms_next_best,delta_rms: ',3f7.3,'  ',
     &        100(a3,'_',a6,' | '))

          write(*,*)'about to do stuff '
          write(*,*)'n1 etc ',n1,nnn,i,j,k,l,nmov(nnn),nmtc(nnn)

          do mmm=1,nmov(nnn)
            xm(mmm)=xlig(nnn,mmm)
            ym(mmm)=ylig(nnn,mmm)
            zm(mmm)=zlig(nnn,mmm)
          enddo

          xf(1)=xxx_mem_region(n1,i)
          yf(1)=yyy_mem_region(n1,i)
          zf(1)=zzz_mem_region(n1,i)
          xf(2)=xxx_mem_region(n1,j)
          yf(2)=yyy_mem_region(n1,j)
          zf(2)=zzz_mem_region(n1,j)
          xf(3)=xxx_mem_region(n1,k)
          yf(3)=yyy_mem_region(n1,k)
          zf(3)=zzz_mem_region(n1,k)
          if(nmtc(nnn).ge.4) then
            xf(4)=xxx_mem_region(n1,l)
            yf(4)=yyy_mem_region(n1,l)
            zf(4)=zzz_mem_region(n1,l)
          endif
          if(nmtc(nnn).ge.5) then
            xf(5)=xxx_mem_region(n1,m)
            yf(5)=yyy_mem_region(n1,m)
            zf(5)=zzz_mem_region(n1,m)
          endif
          if(nmtc(nnn).ge.6) then
            xf(6)=xxx_mem_region(n1,n)
            yf(6)=yyy_mem_region(n1,n)
            zf(6)=zzz_mem_region(n1,n)
          endif
               
C now superimpose

          call pdbsup(nmov(nnn),nmtc(nnn),xm,ym,zm,
     &                xf,yf,zf,rms)
ccc       write(*,*)'RMS_CHECK ',rms

C only write out the non-match atoms...

          do mx=nmtc(nnn)+1,nmov(nnn)
            occ_tmp=real(nligsites)
            beta_tmp=rms_best(nnn)
            nfin_atm=nfin_atm+1
            write(13,851)nfin_atm,string1(nnn,mx)(12:30),
     &                   xm(mx),ym(mx),zm(mx),occ_tmp,beta_tmp
851         format('HETATM',i5,a19,3f8.3,2f6.2)

C also store coords in memory for clash checking later ligands

            xfin(nfin_atm)=xm(mx)
            yfin(nfin_atm)=ym(mx)
            zfin(nfin_atm)=zm(mx)

          enddo

C else we should also write out regions that did NOT fit to any
C template - we might find new motifs this way...

        else

C we first need to find out the best of the shitty rms's we found

          rms_fin=999.9
          nnn_fin=-1
          do nnn=1,nlig
            if(rms_best(nnn).lt.rms_fin) then
              rms_fin=rms_best(nnn)
              nnn_fin=nnn
            endif
          enddo

          write(*,616)initial_pdb(1:30),n1,nnn_fin,rms_fin
616       format('pdbname ',a30,
     &        ' for region# ',i6,' NO GOOD FIT! BEST: ',i6,
     &        ' with rmsd of ',f8.3)

          nfails=nfails+1
          if(nfails.eq.1) then
            open(unit=23,file='ligand_fails.pdb',status='unknown')
          endif

          do mx=1,num_mem_region(n1)
            occ_tmp=real(nfails)
            beta_tmp=rms_best(nnn)
            write(23,591)mx,' SG ',nam_mem_region(n1,mx),
     &                             idr_mem_region(n1,mx),
     &               xxx_mem_region(n1,mx),
     &               yyy_mem_region(n1,mx),
     &               zzz_mem_region(n1,mx),
     &               occ_tmp,beta_tmp
591         format('ATOM ',i6,1x,a4,1x,a3,' A',i4,4x,3f8.3,2f6.2)
          enddo

        endif

C go back and try again if we're not dun yet with this region

        if(iam_not_dun.eq.1) goto 444

C come here if we jumped out earlier with orphans

333     continue

C now process another region

      enddo

C now close the pdb files

      close(13)
      close(23)
      close(31)

C write out our two summary files

      open(unit=21,file=outfile2,status='unknown')
      open(unit=22,file=outfile1,status='unknown')

C figure out how many different types of ligands were identified...

      nunq=0 ! # of types of ligands added
      nadd=0 ! # of ligand molecules added
      do n=1,nlig
        if(ifound(n).eq.1) nunq=nunq+1
        if(nfound(n).ge.1) nadd=nadd+nfound(n)
      enddo
      write(21,'(a30,5x,100i6)')initial_pdb(1:30),(ifound(n),n=1,nlig),
     &                          nunq
      write(22,'(a30,5x,100i6)')initial_pdb(1:30),(nfound(n),n=1,nlig),
     &                          nadd
      close(21)
      close(22)

C write out cysteine statistics

      open(unit=21,file=outfile8,status='unknown')
      num_cys_check=num_cys_singl+num_cys_doubl+
     &              num_cys_lignd+num_cys_orfan
      if(num_cys_check.ne.num_cys_total) then
        write(*,*)
        write(*,*)'WARNING - NON-FATAL ERROR'
        write(*,*)'num_cys_total = ',num_cys_total
        write(*,*)'num_cys_check = ',num_cys_check
        write(*,*)'this can happen when CYS residues are '
        write(*,*)'putative disulfides *AND* are liganded '
        write(*,*)'to metal clusters...'
        write(*,*)
      endif
      write(21,971)initial_pdb(1:30),
     &  num_cys_singl,num_cys_doubl,num_cys_lignd,num_cys_orfan
971   format(a30,5x,' #CYS_singl: ',i6,
     &              ' #CYS_doubl: ',i6,
     &              ' #CYS_lignd: ',i6,
     &              ' #CYS_orfan: ',i6)
      close(21)

C write out a few blank lines to nohup.out

      write(*,*)
      write(*,*)
      write(*,*)

      stop
      end






      subroutine pdbsup(nmov_loc,nmtc_loc,xm_loc,ym_loc,zm_loc,
     &                  xf_loc,yf_loc,zf_loc,rms_loc)

C nmov_loc is the number of total moving atoms (match atoms included)
C nmtc_loc is the number of match atoms

c ----------------------------------------------------------------------PDBS0002
c                               PROGRAM PDBSUP                          PDBS0001
c             Determines rotation matrix and translation vector for     PDBS0002
c              best fit superimposition of two pdb files by solving     PDBS0003
c                      the quaternion eigenvalue problem.               PDBS0004
c                      Code by B.Rupp and S.Parkin (1996)               PDBS0005
c               Lawrence Livermore National Laboratory, BBRP            PDBS0006
c             Method by S.K.Kearsley, Acta Cryst. A45, 208 (1989)       PDBS0007
c              See http://www-structure.llnl.gov for details.           PDBS0008
c ----------------------------------------------------------------------PDBS0009
c !MS$DEBUG                                                             PDBS0010
                                                                        PDBS0011
c --- number of atoms read -                                            PDBS0012

      implicit real(a-h,o-z)
                                                                        PDBS0014
      integer, intent(in)     :: nmov_loc
      integer, intent(in)     :: nmtc_loc
      real, intent(inout)     :: xm_loc(1:nmov_loc)
      real, intent(inout)     :: ym_loc(1:nmov_loc)
      real, intent(inout)     :: zm_loc(1:nmov_loc)
      real, intent(in)        :: xf_loc(1:nmtc_loc)
      real, intent(in)        :: yf_loc(1:nmtc_loc)
      real, intent(in)        :: zf_loc(1:nmtc_loc)
      real,    intent(out)    :: rms_loc

C note that we make rm_loc from xm_loc,ym_loc,zm_loc...

      real                    :: rm_loc(1:nmov_loc,1:3)
      real                    :: rf_loc(1:nmtc_loc,1:3)

      real cack(3)
      real dm(4),vm(4,4),cm(3),cf(3)                                    PDBS0018
      real tr(3),t(3,3),q(4,4)                                          PDBS0019
      real dxp(nmtc_loc,3),dxm(nmtc_loc,3)
      real sumf(3),summ(3)                                               PDBS0021
                                                                        PDBS0022
      deg = 57.29577951                                                 PDBS0023
      zero=10E-30                                                       PDBS0024

      do j=1,nmov_loc
        rm_loc(j,1)=xm_loc(j)
        rm_loc(j,2)=ym_loc(j)
        rm_loc(j,3)=zm_loc(j)
      enddo
      do j=1,nmtc_loc
        rf_loc(j,1)=xf_loc(j)
        rf_loc(j,2)=yf_loc(j)
        rf_loc(j,3)=zf_loc(j)
      enddo
        
      imov=nmtc_loc ! set # atoms for the superpositions... 

c --- initialize all --                                                 PDBS0124
      do i=1,3                                                          PDBS0125
         sumf(i)=0                                                      PDBS0126
         summ(i)=0                                                      PDBS0127
      end do                                                            PDBS0128
      call filmat(4,4,q,0)                                              PDBS0129
                                                                        PDBS0130
c --- sum up all coordinates (in dble precision) to find centre ---     PDBS0131
      do k=1,imov                                                       PDBS0132
         do i=1,3                                                       PDBS0133
            sumf(i)=sumf(i)+rf_loc(k,i)                                     PDBS0134
            summ(i)=summ(i)+rm_loc(k,i)                                     PDBS0135
         end do                                                         PDBS0num_solute_typs
      end do                                                            PDBS0137
      do i=1,3                                                          PDBS0138
         cm(i)=summ(i)/real(imov)                                           PDBS0139
         cf(i)=sumf(i)/real(imov)                                           PDBS0140
         tr(i)=cf(i)-cm(i)                                              PDBS0141
      end do                                                            PDBS0142
                                                                        PDBS0143
c     write(*,'(/a,3f12.3)')' Cen of target molecule  =',(cf(i),i=1,3)  PDBS0144
c     write(*,'(a,3f12.3)') ' Cen of moving molecule  =',(cm(i),i=1,3)  PDBS0145
c     write(*,'(a,3f8.3/)')' T - vector probe -> target =',(tr(i),i=1,3)PDBS0146
                                                                        PDBS0147
c     write (*,'(a)')' Creating coordinate differences.......'          PDBS0148
c --- create coordinate differences delta x plus (dxp) and minus (dxm)  PDBS0149
      do k=1,imov                                                       PDBS0150
         do j=1,3                                                       PDBS0151
            dxm(k,j)=rm_loc(k,j)-cm(j)-(rf_loc(k,j)-cf(j))                      PDBS0152
            dxp(k,j)=rm_loc(k,j)-cm(j)+(rf_loc(k,j)-cf(j))                      PDBS0153
         end do                                                         PDBS0154
      end do                                                            PDBS0155
                                                                        PDBS0156
c --- fill upper triangle of (symmetric) quaternion matrix --           PDBS0157
c     write (*,'(a)')' Filling quaternion matrix ............'          PDBS0158
      do k=1,imov                                                       PDBS0159
c ---    diags are sums of squared cyclic coordinate differences        PDBS0160
         q(1,1)=q(1,1)+dxm(k,1)**2+dxm(k,2)**2+dxm(k,3)**2              PDBS0161
         q(2,2)=q(2,2)+dxp(k,2)**2+dxp(k,3)**2+dxm(k,1)**2              PDBS0162
         q(3,3)=q(3,3)+dxp(k,1)**2+dxp(k,3)**2+dxm(k,2)**2              PDBS0163
         q(4,4)=q(4,4)+dxp(k,1)**2+dxp(k,2)**2+dxm(k,3)**2              PDBS0164
c ---    cross differences                                              PDBS0165
         q(1,2)=q(1,2)+dxp(k,2)*dxm(k,3)-dxm(k,2)*dxp(k,3)              PDBS0166
         q(1,3)=q(1,3)+dxm(k,1)*dxp(k,3)-dxp(k,1)*dxm(k,3)              PDBS0167
         q(1,4)=q(1,4)+dxp(k,1)*dxm(k,2)-dxm(k,1)*dxp(k,2)              PDBS0168
         q(2,3)=q(2,3)+dxm(k,1)*dxm(k,2)-dxp(k,1)*dxp(k,2)              PDBS0169
         q(2,4)=q(2,4)+dxm(k,1)*dxm(k,3)-dxp(k,1)*dxp(k,3)              PDBS0170
         q(3,4)=q(3,4)+dxm(k,2)*dxm(k,3)-dxp(k,2)*dxp(k,3)              PDBS0171
      end do                                                            PDBS0172
c --- fill the rest by transposing it onto itself                       PDBS0173
      call trpmat(4,q,q)                                                PDBS0174
c     write (*,'(/a)')                                                  PDBS0175
c    &     '       q(1)         q(2)         q(3)        q(4)'          PDBS0176
c     do i=1,4                                                          PDBS0177
c        write(*,'(4e13.5)') (q(i,j),j=1,4)                             PDBS0178
c     end do                                                            PDBS0179
                                                                        PDBS0180
c --- orthogonalization by jacobi rotation = solution of EV -problem -- PDBS0181
c     write (*,'(/a)')' Jacobi orthogonalization ..........'            PDBS0182
      n=4                                                               PDBS0183
      ns=4                                                              PDBS0184
      call jacobi(q,n,ns,dm,vm,nmrot)                                   PDBS0185
c --- sort eigenvectors after eigenvalues, descending --                PDBS0186
c     write (*,'(a/)')' Sorting eigenvalues/vectors .......'            PDBS0187
      call eigsrt(dm,vm,n,ns)                                           PDBS0188
c     write (*,'(a,i2,a)')' Eigenvalues and Eigenvectors (',            PDBS0189
c    & nmrot,' Jacobi rotations)'                                       PDBS0190
c     write (*,'(a)') '      e(1)        e(2)        e(4)        e(4)'  PDBS0191
c     write (*,'(4e12.5,i5)') (dm(j),j=1,4)                             PDBS0192
c     write (*,'(a)') '      ev(1)       ev(2)       ev(3)       ev(4)' PDBS0193
c     do i=1,4                                                          PDBS0194
c        write(*,'(4f12.6)') (vm(i,j),j=1,4)                            PDBS0195
c     end do                                                            PDBS0196
                                                                        PDBS0197
c --- the smallest eigenvector contains best fit srs                    PDBS0198
      rmsd=sqrt(abs(dm(4)/imov))                                        PDBS0199

c     write(*,*)'my rmsd = ',rmsd
c     write(*,'(/a/)')                                                  PDBS0200
c    & ' The smallest eigenvalue represents s.r.s. of best fit'         PDBS0201
c     write(*,'(a)')                                                    PDBS0202
c    & ' Constructing the best fit rotation matrix from associated'     PDBS0203
c     write(*,'(a/)') ' eigenvector elements (last column).....'        PDBS0204
                                                                        PDBS0205
c --- fill the rotation matrix which is made of elements from 4th EV    PDBS0206
      t(1,1)=vm(1,4)**2+vm(2,4)**2-vm(3,4)**2-vm(4,4)**2                PDBS0207
      t(2,1)=2*(vm(2,4)*vm(3,4)+vm(1,4)*vm(4,4))                        PDBS0208
      t(3,1)=2*(vm(2,4)*vm(4,4)-vm(1,4)*vm(3,4))                        PDBS0209
      t(1,2)=2*(vm(2,4)*vm(3,4)-vm(1,4)*vm(4,4))                        PDBS0210
      t(2,2)=vm(1,4)**2+vm(3,4)**2-vm(2,4)**2-vm(4,4)**2                PDBS0211
      t(3,2)=2*(vm(3,4)*vm(4,4)+vm(1,4)*vm(2,4))                        PDBS0212
      t(1,3)=2*(vm(2,4)*vm(4,4)+vm(1,4)*vm(3,4))                        PDBS0213
      t(2,3)=2*(vm(3,4)*vm(4,4)-vm(1,4)*vm(2,4))                        PDBS0214
      t(3,3)=vm(1,4)**2+vm(4,4)**2-vm(2,4)**2-vm(3,4)**2                PDBS0215
                                                                        PDBS0216
c     do i=1,3                                                          PDBS0217
c        write(*,'(3f11.5)') (t(i,j),j=1,3)                             PDBS0218
c     end do                                                            PDBS0219
                                                                        PDBS0220
c --- reset dxm to store the individual rmsd's in it now -              PDBS0221
      call filmat(nmtc_loc,3,dxm,0)                                          PDBS0222
                                                                        PDBS0223
c --- xm and xf are not translated                                      PDBS0224
      do k=1,nmov_loc                                                    PDBS0225
c ---    subtract cm                                                    PDBS0226
         rm_loc(k,1)=rm_loc(k,1)-cm(1)                                             PDBS0227
         rm_loc(k,2)=rm_loc(k,2)-cm(2)                                             PDBS0227
         rm_loc(k,3)=rm_loc(k,3)-cm(3)                                             PDBS0227
C
C AHE
C make call to rotvec easier
C
         cack(1)=rm_loc(k,1)
         cack(2)=rm_loc(k,2)
         cack(3)=rm_loc(k,3)
c ---    rotate it                                                      PDBS0228
         call rotvec(3,cack,t)                                          PDBS0229
c ---    now add cf                                                     PDBS0230
         xm_loc(k)=cack(1)
         ym_loc(k)=cack(2)
         zm_loc(k)=cack(3)
         xm_loc(k)=xm_loc(k)+cf(1)                                          PDBS0231
         ym_loc(k)=ym_loc(k)+cf(2)                                          PDBS0231
         zm_loc(k)=zm_loc(k)+cf(3)                                          PDBS0231
c        do i=1,3                                                       PDBS0232
c           dxm(k,i)=sqrt((rf_loc(k,i)-rm_loc(k,i))**2)                         PDBS0233
c        end do                                                         PDBS0234
      end do                                                            PDBS0235
                                                                        PDBS0236
      rms_loc=rmsd
c     write(6,3)s                                                       PDBS0260
                                                                        PDBS0293
 0003      format(a)

      end                                                               PDBS0354
                                                                        PDBS0355
      subroutine inkey (answ)                                           INKE0001
c ----------------------------------------------------------------------INKE0002
c     reads a key as answer                                             INKE0003
c ----------------------------------------------------------------------INKE0004
      character answ                                                    INKE0005
      read (*,'(a1)') answ                                              INKE0006
      call upstrg (answ,1)                                              INKE0007
      if ((answ.eq.' ').or.(answ.eq.char(13))) answ='Y'                 INKE0008
      return                                                            INKE0009
      end                                                               INKE0010
                                                                        INKE0011
      subroutine upstrg(strg,istrln)                                    UPST0001
c ----------------------------------------------------------------------UPST0002
c converts string str$ of lenght istrlen to upcase                      UPST0003
c ----------------------------------------------------------------------UPST0004
      character strg(istrln)                                            UPST0005
      integer      iascii                                               UPST0006
C-----change to upper case                                              UPST0007
         do 3031 i=1,istrln                                             UPST0008
            if (ichar(strg(i)).ge.95.and.ichar(strg(i)).le.122) then    UPST0009
               iascii=ichar(strg(i))-32                                 UPST0010
               strg(i)=char(iascii)                                     UPST0011
            endif                                                       UPST0012
 3031    continue                                                       UPST0013
      return                                                            UPST0014
      end                                                               UPST0015
                                                                        UPST0016
      subroutine trpmat(n,t,tr)                                         TRPM0001
c --- transpose matrix -------------------------------------------------TRPM0002
      real t(n,n), tr(n,n)                                              TRPM0003
      do i=1,n                                                          TRPM0004
         do j=1,n                                                       TRPM0005
            tr(j,i)=t(i,j)                                              TRPM0006
         end do                                                         TRPM0007
      end do                                                            TRPM0008
      return                                                            TRPM0009
      end                                                               TRPM0010
                                                                        TRPM0011
      subroutine filmat(n,m,r,ifil)                                     FILM0001
      real r(n,m)                                                       FILM0002
c --- initialize matrix ------------------------------------------------FILM0003
      do i=1,n                                                          FILM0004
         do j=1,m                                                       FILM0005
            r(i,j)=ifil                                                 FILM0006
         end do                                                         FILM0007
      end do                                                            FILM0008
      return                                                            FILM0009
      end                                                               FILM0010
                                                                        FILM0011
      subroutine rotvec (n,v,t)                                         ROTV0001
c --- multiply vector with matrix --------------------------------------ROTV0002
      real t(n,n), v(n),s(n)                                            ROTV0003
                                                                        ROTV0004
      do i=1,n                                                          ROTV0005
         s(i)=v(i)                                                      ROTV0006
         v(i)=0.0                                                       ROTV0007
      end do                                                            ROTV0008
      do i=1,n                                                          ROTV0009
         do j=1,n                                                       ROTV0010
            v(i)=v(i)+s(j)*t(i,j)                                       ROTV0011
         end do                                                         ROTV0012
      end do                                                            ROTV0013
      return                                                            ROTV0014
      end                                                               ROTV0015
                                                                        ROTV0016
      SUBROUTINE eigsrt(d,v,n,np)                                       EIGS0001
c ----------------------------------------------------------------------EIGS0002
      INTEGER n,np                                                      EIGS0003
      REAL d(np),v(np,np)                                               EIGS0004
      INTEGER i,j,k                                                     EIGS0005
      REAL p                                                            EIGS0006
      do 13 i=1,n-1                                                     EIGS0007
        k=i                                                             EIGS0008
        p=d(i)                                                          EIGS0009
        do 11 j=i+1,n                                                   EIGS0010
          if(d(j).ge.p)then                                             EIGS0011
            k=j                                                         EIGS0012
            p=d(j)                                                      EIGS0013
          endif                                                         EIGS0014
11      continue                                                        EIGS0015
        if(k.ne.i)then                                                  EIGS0016
          d(k)=d(i)                                                     EIGS0017
          d(i)=p                                                        EIGS0018
          do 12 j=1,n                                                   EIGS0019
            p=v(j,i)                                                    EIGS0020
            v(j,i)=v(j,k)                                               EIGS0021
            v(j,k)=p                                                    EIGS0022
12        continue                                                      EIGS0023
        endif                                                           EIGS0024
13    continue                                                          EIGS0025
      return                                                            EIGS0026
      END                                                               EIGS0027
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             EIGS0028
                                                                        EIGS0029
      SUBROUTINE jacobi(a,n,np,d,v,nrot)                                JACO0001
c ----------------------------------------------------------------------JACO0002
c     modified from numerical recipes book                              JACO0003
c     one needs to set the threshold for sm from sm.eq.0 to sm.lt.10E-30JACO0004
c     (anything in this range would be ok) due to underflow errors on   JACO0005
c     some computers/compilers.                                         JACO0006
c ----------------------------------------------------------------------JACO0007
      PARAMETER (nmax=500)                                              JACO0008
                                                                        JACO0009
      INTEGER n,np,nrot                                                 JACO0010
      REAL a(np,np),d(np),v(np,np)                                      JACO0011
      INTEGER i,ip,iq,j,maxrot                                          JACO0012
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),zero          JACO0013
                                                                        JACO0014
c --- zero set and iteration maximum                                    JACO0015
      zero=10E-30                                                       JACO0016
      maxrot=50                                                         JACO0017
                                                                        JACO0018
      do 12 ip=1,n                                                      JACO0019
        do 11 iq=1,n                                                    JACO0020
          v(ip,iq)=0.                                                   JACO0021
11      continue                                                        JACO0022
        v(ip,ip)=1.                                                     JACO0023
12    continue                                                          JACO0024
      do 13 ip=1,n                                                      JACO0025
        b(ip)=a(ip,ip)                                                  JACO0026
        d(ip)=b(ip)                                                     JACO0027
        z(ip)=0.                                                        JACO0028
13    continue                                                          JACO0029
      nrot=0                                                            JACO0030
      do 24 i=1,maxrot                                                  JACO0031
        sm=0.                                                           JACO0032
        do 15 ip=1,n-1                                                  JACO0033
          do 14 iq=ip+1,n                                               JACO0034
            sm=sm+abs(a(ip,iq))                                         JACO0035
14        continue                                                      JACO0036
15      continue                                                        JACO0037
c ---   modified convergence threshold ---                              JACO0038
        if(sm.lt.zero)return                                            JACO0039
        if(i.lt.4)then                                                  JACO0040
          tresh=0.2*sm/n**2                                             JACO0041
        else                                                            JACO0042
          tresh=0.                                                      JACO0043
        endif                                                           JACO0044
        do 22 ip=1,n-1                                                  JACO0045
          do 21 iq=ip+1,n                                               JACO0046
            g=100.*abs(a(ip,iq))                                        JACO0047
            if((i.gt.4).and.(abs(d(ip))+                                JACO0048
     &g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then            JACO0049
              a(ip,iq)=0.                                               JACO0050
            else if(abs(a(ip,iq)).gt.tresh)then                         JACO0051
              h=d(iq)-d(ip)                                             JACO0052
              if(abs(h)+g.eq.abs(h))then                                JACO0053
                t=a(ip,iq)/h                                            JACO0054
              else                                                      JACO0055
                theta=0.5*h/a(ip,iq)                                    JACO0056
                t=1./(abs(theta)+sqrt(1.+theta**2))                     JACO0057
                if(theta.lt.0.)t=-t                                     JACO0058
              endif                                                     JACO0059
              c=1./sqrt(1+t**2)                                         JACO0060
              s=t*c                                                     JACO0061
              tau=s/(1.+c)                                              JACO0062
              h=t*a(ip,iq)                                              JACO0063
              z(ip)=z(ip)-h                                             JACO0064
              z(iq)=z(iq)+h                                             JACO0065
              d(ip)=d(ip)-h                                             JACO0066
              d(iq)=d(iq)+h                                             JACO0067
              a(ip,iq)=0.                                               JACO0068
              do 16 j=1,ip-1                                            JACO0069
                g=a(j,ip)                                               JACO0070
                h=a(j,iq)                                               JACO0071
                a(j,ip)=g-s*(h+g*tau)                                   JACO0072
                a(j,iq)=h+s*(g-h*tau)                                   JACO0073
16            continue                                                  JACO0074
              do 17 j=ip+1,iq-1                                         JACO0075
                g=a(ip,j)                                               JACO0076
                h=a(j,iq)                                               JACO0077
                a(ip,j)=g-s*(h+g*tau)                                   JACO0078
                a(j,iq)=h+s*(g-h*tau)                                   JACO0079
17            continue                                                  JACO0080
              do 18 j=iq+1,n                                            JACO0081
                g=a(ip,j)                                               JACO0082
                h=a(iq,j)                                               JACO0083
                a(ip,j)=g-s*(h+g*tau)                                   JACO0084
                a(iq,j)=h+s*(g-h*tau)                                   JACO0085
18            continue                                                  JACO0086
              do 19 j=1,n                                               JACO0087
                g=v(j,ip)                                               JACO0088
                h=v(j,iq)                                               JACO0089
                v(j,ip)=g-s*(h+g*tau)                                   JACO0090
                v(j,iq)=h+s*(g-h*tau)                                   JACO0091
19            continue                                                  JACO0092
              nrot=nrot+1                                               JACO0093
            endif                                                       JACO0094
21        continue                                                      JACO0095
22      continue                                                        JACO0096
        do 23 ip=1,n                                                    JACO0097
          b(ip)=b(ip)+z(ip)                                             JACO0098
          d(ip)=b(ip)                                                   JACO0099
          z(ip)=0.                                                      JACO0100
23      continue                                                        JACO0101
24    continue                                                          JACO0102
ccc   pause 'too many iterations in jacobi'                             JACO0103
      return                                                            JACO0104
      END                                                               JACO0105
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             JACO0106



         SUBROUTINE init_random_seed()

            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)

            CALL RANDOM_SEED(PUT = seed)
            DEALLOCATE(seed)

         END SUBROUTINE


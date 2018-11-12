;.....GENERAZIONE DEI PROFILI DI TUTTO IL CATALOGO A TUTTE LE SCALE DI FILTRAGGIO ....

pro do_all

tini=systime(1)
 fold_out='/home/STUDENTI/matteo.angelinelli/Output/'
 fold='/home/data/DATA/ISC/'
  cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]   ;...list of available clusters (considering z=0 snapshots)

  sn_list=[199,196,195,193,199,243,228,242,102,100,115,167,106,109] ;....snapshot number for z=0 (actually, last available snapshot, very close to z=0)
                                ;....clusters IT1,IT3,IT7,IT10 belongs
                                ;to a first generation of the ITASCA
                                ;sample, with a coarser time-sampling
                                ;for the outputs (hence their snapshot
                                ;number for z=0 is smaller)

  nc=size(cl_list)
  print,"going to process ",nc(1),"clusters"

deltam=fltarr(nc(1))
alpha=fltarr(nc(1))

close,3
close,5
;openw,5,fold_out+'catalog_cluster_OK.dat'
;printf,5,'Cluster -- R200 [Mpc] -- M200_hyd [Mo] -- M200_tot [Mo] -- DeltaM/M200_tot -- Alpha200',format='(A)'

alfa1=0.
deltam1=0.
  
 for c=0,nc(1)-1 do begin
;c=3
   nt_pressure,fold,cl_list(c),sn_list(c),alfa1,deltam1    ;..we just call the analysis routine below for each cluster

   ;nt_pressure_R200,fold,cl_list(c),sn_list(c)

deltam(c)=deltam1
alpha(c)=alfa1

 endfor   

close,4
openw,4,fold_out+'alpha_parameter_OK.dat'
printf,4,'Cluster -- Alpha200 -- (Mtot-Mhyd)/Mtot',format='(A)'
for c=0,nc(1)-1 do begin
printf,4,cl_list(c),alpha(c),deltam(c),format='(A,F,F)'
endfor
close,5
close,4

;do_norm
;do_norm_r200_3   

print, '...procedure complete in ',systime(1)-tini, ' seconds...'
end

;...GENERAZIONE DEL CATALOGO FILTRATTO A SCALE FISSE.....

pro NT_pressure,foldc,cluster,snapc,alfa,deltamm

  t0=systime(1)

  fold_out='/home/STUDENTI/matteo.angelinelli/Output/'    ;...output folder

  CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine
 
  fold=foldc+cluster+'/'
  snap='0'+string(snapc,'(i3)')
  file1=fold+'/declust_z'+snap   ;...the field names are Density, Temperature, Dark_Matter_Densiy
  file2=fold+'declust_v_z'+snap ;..." x-velocity, y-velocity, z-velocity

  if cluster eq 'IT3' or cluster eq 'IT7' or cluster eq 'IT6' or cluster eq 'IT15B' then file2=file1  ;...these files where written in a slightly different way

  fileconv='deDD'+snap+'.conv2' ;...conversion factors from enzo's internal units to cgs

  openr,3,fold+fileconv
  readf,3,snapshot_id
  readf,3,time
  readf,3,redshift
  readf,3,convd  ;...conversion factor from enzo's to g/cm^3
  readf,3,convv   ;...      "      "        cm/s
  close,3

  print,cluster,redshift,convd,convv

  ;....the conversion factors change with redshift and simulation, must be found into enzo output files or are written elsewhere during the data reconstruction (ask Franco)
 ;....HDF5 (see if it's easy to download the hdf5 libraries for shell commands
  ;...https://www.hdfgroup.org/downloads/


  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Density')
  datasp=h5d_get_space(itemp)
  dens= H5D_READ(itemp)
  h5d_close,itemp

  
  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Dark_Matter_Density')
  datasp=h5d_get_space(itemp)
  dens_dm = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'Temperature')
  datasp=h5d_get_space(itemp)
  temp = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido


  ido=h5f_open(file2)
  itemp=h5d_open(ido,'x-velocity')
  datasp=h5d_get_space(itemp)
  vx = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'y-velocity')
  datasp=h5d_get_space(itemp)
  vy = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'z-velocity')
  datasp=h5d_get_space(itemp)
  vz = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido

;....some manipulation of grid size needed only in the rare case the
;HDF5 dataset is not a cube

  n3=size(dens)                    ;...this gives the size of the density field just read
  print,n3

  n0=n3(1)
  n1=n3(2)
  n2=n3(3)
   n=n0
  if n0 ne n1 or n0 ne n2 or n1 ne n2 then begin
     print,"making the dataset cubic"
                                ;....probably there is a much more
                                ;elegant way of doing the following, anyways....
     n=min(n3(1:3))
     dnew=dens(0:n-1,0:n-1,0:n-1)
     tnew=temp(0:n-1,0:n-1,0:n-1)
     dmnew=dens_dm(0:n-1,0:n-1,0:n-1)
     vxnew=vx(0:n-1,0:n-1,0:n-1)
     vynew=vy(0:n-1,0:n-1,0:n-1)
     vznew=vz(0:n-1,0:n-1,0:n-1)
     dens=dnew
     dens_dm=dmnew
     temp=tnew
     vx=vxnew
     vy=vynew
     vz=vznew
     vxnew=0
     vynew=0
     vznew=0
  end   

if cluster eq 'IT10' then begin 
res=10 ;...kpc
endif
if cluster ne 'IT10' then begin
res=20
endif

dens=dens*convd
dens_dm=dens_dm*convd

t1=systime(1)
print, '...read complete in',(systime(1)-t0),' seconds...'

vol=3.*alog10(res)+3.*alog10(3.086*10.^21.)
;...the first 3 are the input un-filtered velocities, assumed to be in cgs
;...the input density must be in cgs
;...vol=log10(cell volume in cgs)
;...resdhift
;....flux is the (output) kinetic energy flux through shocked cells, in log10(erg/s)
;....macx is the (output) Mach number of shocked cells
;....div is the (output) 3D divergence of the velocity field, in cm/s/cell

;a=shock_finder(convv*vx,convv*vy,convv*vz,dens,temp,vol-3.*alog10(1.+redshift),flux,macx,div) 

t2=systime(1)
print, '...shocks finder complete in',(systime(1)-t1),' seconds...'
  
  xm=where(dens+dens_dm eq max(dens+dens_dm))
  xc=array_indices(dens,xm)

  print,"the center of this cluster is",xc

  rx=indgen(n)-xc(0)
  ry=indgen(n)-xc(1)
  rz=indgen(n)-xc(2)
  rx=reform(rx,n,1,1)
  ry=reform(ry,1,n,1)
  rz=reform(rz,1,1,n)
  ra=rebin(rx,n,n,n)^2
  ra=temporary(ra)+rebin(ry,n,n,n)^2                 ; conserve memory
  ra=temporary(ra)+rebin(rz,n,n,n)^2
  ra=fix(temporary(sqrt(ra)))  ;...array of 3D radial distances


print, '...grid complete in ',systime(1)-t2,' seconds...'

n=320.
cont_scale=0
scale=[3,5,10,20,30]
return_scale:
t3=systime(1)
sc=scale(cont_scale)
print, sc

  vx_smooth=smooth(vx,sc)
  vy_smooth=smooth(vy,sc)
  vz_smooth=smooth(vz,sc)

  vturbx=vx-vx_smooth
  vturby=vy-vy_smooth
  vturbz=vz-vz_smooth

  vmod_turb=convv*sqrt(vturbx^2.+vturby^2.+vturbz^2.)  ;...velocity module in cm/s

vx_smooth=0
vy_smooth=0
vz_smooth=0
vturbx=0
vturby=0
vturbz=0

t4=systime(1)
print, '...velocity fields complete in',systime(1)-t3,' seconds...'


pd2=fltarr(max(ra))
pdm2=fltarr(max(ra))
pd=fltarr(max(ra))
pdm=fltarr(max(ra))
ptemp=fltarr(max(ra))
pgas=fltarr(max(ra))
pmass=fltarr(max(ra))
pnt=fltarr(max(ra)) 
ptot=fltarr(max(ra))
 ratio=fltarr(max(ra))
  kb=1.38*10.^(-16.)
  mh=1.67*10.^(-24.)
  mu=1.1
  mue=0.59

  for r=0,max(ra)-1 do begin
  
   rw=where(ra ge r and ra lt r+1,nw)   ;…we can also do r
   dw=percentiles(dens(rw),value=[0.9])


 ; rk=where(ra ge r and ra lt r+1 and macx lt 1.3 and dens lt dw,nr)   ;…we can also do r
; ry=where(ra ge r and ra lt r+1 and macx lt 1.3,ny)   ;…we can also do r 
  rz=where(ra ge r and ra lt r+1 and dens lt dw,nz)

;print,(float(nr)/float(ny))*100.,(float(nz)/float(ny))*100.,(float(nr)/float(nw))*100.,(float(nz)/float(nw))*100.,(float(ny)/float(nw))*100.
;print,(float(nz)/float(nw))*100.
  ;  pnt(r)=mean((1./3.)*dens(rk)*vmod_turb(rk)^(2.))  
  ; pgas(r)=mean(kb*temp(rk)*dens(rk)/float(mh*mue))
  ;  ptot(r)=pgas(r)+pnt(r)
  ;  ratio(r)=pnt(r)/ptot(r)    

    pmass(r)=mean(kb*temp(rz)*dens(rz)/float(mh*mue))
    pd(r)=mean(dens(rw))
    pdm(r)=mean(dens_dm(rw))
    pd2(r)=mean(dens(rz))
    pdm2(r)=mean(dens_dm(rz))
    ptemp(r)=mean(temp(rz))
    
endfor

vol=3.*alog10(res)+3.*alog10(3.086*10.^21.)
pgas=alog10(pgas)+vol
pnt=alog10(pnt)+vol
pmass=alog10(pmass)+vol

x=res*(0.5+indgen(max(ra)))

overd=200 ;....overdensity 
r200=rvir(x,pd,pdm,overd)   ;...procedure to derive the radius enclosing a >overd overdensity
x=x/float(r200*res)
rhoc=10.^(-29.063736)
pd=pd*rhoc
pdm=pdm*rhoc

print, '...profile complete in',systime(1)-t4,' seconds...'

sc_name=['3','5','10','20','30']  

;set_plot,'x'
;window,0
;plot,x,pd,/xlog,/ylog
;window,1
;plot,x,pdm,/xlog,/ylog
;window,2
;plot,x,ptemp,/xlog,/ylog
;window,3
;plot,x,ratio,/ylog,/xlog
;window,4
;plot,x,pnt,/xlog,yrange=[45.,60.]
;window,5
;plot,x,pgas,/xlog,yrange=[45.,60.]
;window,6 
;plot,x,pmass,/xlog,yrange=[45.,60.]

;  close,1
; openw,1,fold_out+'Prof/'+cluster+'_pressure_'+sc_name(cont_scale)+'_profile.dat'
;printf,1,max(ra)
 ; for i=0,max(ra)-1 do begin
 ; printf,1,x(i),pnt(i),pgas(i),ratio(i),pd(i),pdm(i),format='(E,E,E,E,E,E)'
 ; endfor
 ; close,1

cont_scale=cont_scale+1

if (cont_scale eq 1) then begin
;openw,2,fold_out+'profile_d_dm_t_'+cluster+'_OK.dat'
;printf,2,'R/R200 -- Density [g cm^-3] -- DM Density [g cm^-3] -- Temperature [K]' 
;for i =0,max(ra)-1 do begin
; printf,2,x(i),pd(i),pdm(i),ptemp(i),format='(E,E,E,E)'
;endfor
;close,2

    kpc_cm=3.086*10.^(21.)
    G=6.67*1e-8
    pmass=10.^(pmass-vol)
    msol=1.9e33
    p200=pmass(r200)*1e30
    rho200=pd2(r200)*1e30

    m_hyd=mass_hydP(x,pmass,pd2,r200*res,res)
 
    m_tot=mass_tot(res,r200,overd)

    dadr=deriv(ratio)
    dadr/=(2.*res*kpc_cm)
    dadr200=dadr(r200)
    deltam=10.^(m_tot-40)-10.^(m_hyd-40.)
    mdm=10.^(m_tot-40.)+deltam 

   alfa=1.-((rho200*G*10.^(m_hyd-40.)+sqrt(abs(rho200*G*10.^(m_hyd-40.)-4.*rho200*G*10.^(m_tot-40.)*p200*res^2.*r200^2.*(kpc_cm*1e-20)^2.*dadr200)))/(2.*10.^(m_tot-40.)*rho200*G))

   ;deltamm=deltam/10.^(m_tot-40.)
    deltamm=1.-10.^(m_hyd-m_tot)

print,10.^(m_hyd-alog10(msol)),10.^(m_tot-alog10(msol)),alfa,deltamm,r200*res

; printf,5,cluster+'   ',r200*res/1.e3,10.^(m_hyd-alog10(msol)),10.^(m_tot-alog10(msol)),deltamm,alfa,format='(A,F,E,E,E,E)'
endif
;if (cont_scale le 4) then begin 
;goto,return_scale
;endif
;if (cont_scale eq 5) then begin 

;t3=systime(1)

;sc=r200/3.

;print, sc

;vtot=fltarr(n,n,n)
;vturb=fltarr(n,n,n)

;  vx_smooth=smooth(vx,sc)
;  vy_smooth=smooth(vy,sc)
 ; vz_smooth=smooth(vz,sc)

 ; vturbx=vx-vx_smooth
 ; vturby=vy-vy_smooth
 ; vturbz=vz-vz_smooth

 ; vmod_turb=convv*sqrt(vturbx^2.+vturby^2.+vturbz^2.)  ;...velocity module in cm/s

vx_smooth=0
vy_smooth=0
vz_smooth=0
vturbx=0
vturby=0
vturbz=0

;print, '...velocity fields complete in',systime(1)-t3,' seconds...'

;t4=systime(1)

;pgas=fltarr(max(ra))
;pnt=fltarr(max(ra)) 
;ptot=fltarr(max(ra))
;ratio=fltarr(max(ra))

  ;for r=0,max(ra)-1 do begin

 ; rw=where(ra ge r and ra lt r+1,nw)   ;…we can also do r
 ; dw=percentiles(dens(rw),value=[0.9]) 
 ; rk=where(ra ge r and ra lt r+1 and macx lt 1.3 and dens lt dw,nr)   ;…we can also do r
 ; rz=where(ra ge r and ra lt r+1 and dens lt dw,nz)

 ;   pnt(r)=mean((1./3.)*dens(rk)*vmod_turb(rk)^(2.))  
 ;   pgas(r)=mean(kb*temp(rk)*dens(rk)/float(mh*mue))
 ;   ptot(r)=pgas(r)+pnt(r)
 ;   pd(r)=mean(dens(rz))
 ;   pdm(r)=mean(dens_dm(rz))
 ;   ratio(r)=pnt(r)/ptot(r)
;endfor


;vol=3.*alog10(res)+3.*alog10(3.086*10.^21.)
;pgas=alog10(pgas)+vol
;pnt=alog10(pnt)+vol

;x=res*(0.5+indgen(max(ra)))

;overd=200 ;....overdensity 
;r200=rvir(x,pd,pdm,overd)   ;...procedure to derive the radius enclosing a >overd overdensity
;x=x/float(r200*res)
;rhoc=10.^(-29.063736)
;pd=pd*rhoc
;pdm=pdm*rhoc

print, '...profile complete in',systime(1)-t4,' seconds...'
   
;  close,6
;  openw,6,fold_out+'Prof/'+cluster+'_'+'scalamobile_pressure_profile.dat'
; for i=0,max(ra)-1 do begin
;  printf,6,x(i),pnt(i),pgas(i),ratio(i),pd(i),pdm(i),format='(E,E,E,E,E,E)'
;  endfor
;  close,6

vx=0
vy=0
vz=0 
dens=0
temp=0
dens_dm=0

;endif

end


;.....CALCOLO DEL RAGGIO DI OVERDENSITY.....


function rvir,r,d,dm,overd,mass
rhoc=10.^(-29.063736)    ;….<———aggiunto
d/=rhoc                  ;..<———aggiunto
dm/=rhoc                 ;..<———aggiunto
nr=size(r)
vshell=4.*!pi*indgen(1+nr(1))^2. ;...volume of each shell in cells
dens_old=1e6
for i=0,nr(1)-1 do begin
dens1=total(vshell(0:i)*(d(0:i)+dm(0:i)))/float((total(vshell(0:i)))) ;...enclosed mean density (units are already in critical density units)
if dens1 lt overd and dens_old gt overd then goto,endr
dens_old=dens1
endfor

endr:
radius=i


return,radius
end

;.....SHOCK FINDER DA VAZZA&JONES 2017.....

function shock_finder,velx,vely,velz,dens,temp,vcell,flux,macx,div


;...for details look into Vazza, Jones+ 2017 MNRAS, Sec.3.1
  mthr=1.3 ;.....low Mach number limiter 

  electronmass=9.109E-27
  protonmass=1.6726E-24  ;in cgs
  kboltz=1.381E-16  ;...in cgs
  molweight=1
  angamma=1.66667
  mw_pm=1/float(protonmass)

  n3=size(dens)
  n=n3(1)

  div=fltarr(n,n,n)
  macx=fltarr(n,n,n)
  macy=fltarr(n,n,n)
  macz=fltarr(n,n,n)

  nsh=uint(0)

  for c=1,n-2 do begin
    for b=1,n-2 do begin
      for a=1,n-2 do begin
        div(a,b,c)=0.5*(velx(a+1,b,c)-velx(a-1,b,c)+vely(a,b+1,c)-vely(a,b-1,c)+velz(a,b,c+1)-velz(a,b,c-1))
      endfor
    endfor
  endfor

  v2=fltarr(2)

  fac=angamma*kboltz*mw_pm

  i3=where(div lt 0 and temp gt 1e5 ,nn)
  print,nn,'candidate shocks'
  bo=array_indices(div,i3)
  print,bo(*,0)

  for i=0L,nn-1 do begin
    a=bo(0,i)
    b=bo(1,i)
    c=bo(2,i)

    dvx=-1*(velx(a-1,b,c)-velx(a+1,b,c))  ;...we need velocity in cm/s here --> k\m/s
    dvy=-1*(vely(a,b-1,c)-vely(a,b+1,c))
    dvz=-1*(velz(a,b,c-1)-velz(a,b,c+1))

    if dvx lt 0 and temp(a+1,b,c) gt temp(a-1,b,c) then begin
      mx=abs(dvx)/float(sqrt(fac*temp(a-1,b,c)))
      mx=(4*mx+sqrt(16*mx*mx+36))*0.166666
      v2(0)=mx
      v2(1)=macx(a+1,b,c)
      macx(a+1,b,c)=max(v2)
    endif else begin
      if dvx lt 0 and temp(a+1,b,c) lt temp(a-1,b,c) then begin
        mx=abs(dvx)/float(sqrt(fac*temp(a+1,b,c)))
        mx=(4*mx+sqrt(16*mx*mx+36))*0.166666
        v2(0)=mx
        v2(1)=macx(a-1,b,c)
        macx(a-1,b,c)=max(v2)
      endif
    endelse


    if dvy lt 0 and temp(a,b+1,c) gt temp(a,b-1,c) then begin
      my=abs(dvy)/float(sqrt(fac*temp(a,b-1,c)))
      my=(4*my+sqrt(16*my*my+36))*0.166666
      v2(0)=my
      v2(1)=macy(a,b+1,c)
      macy(a,b+1,c)=max(v2)
    endif else begin
      if dvy lt 0 and temp(a,b+1,c) lt temp(a,b-1,c) then begin
        my=abs(dvy)/float(sqrt(fac*temp(a,b+1,c)))
        my=(4*my+sqrt(16*my*my+36))*0.166666
        v2(0)=my
        v2(1)=macy(a,b-1,c)
        macy(a,b-1,c)=max(v2)
      endif
    endelse

    if dvz lt 0 and temp(a,b,c+1) gt temp(a,b,c-1) then begin
      mz=abs(dvz)/float(sqrt(fac*temp(a,b,c-1)))
      mz=(4*mz+sqrt(16*mz*mz+36))*0.166666
      v2(0)=mz
      v2(1)=macz(a,b,c+1)
      macz(a,b,c+1)=max(v2)
    endif else begin
      if dvz lt 0 and temp(a,b,c+1) lt temp(a,b,c-1) then begin
        mz=abs(dvz)/float(sqrt(fac*temp(a,b,c+1)))
        mz=(4*mz+sqrt(16*mz*mz+36))*0.166666
        v2(0)=mz
        v2(1)=macz(a,b-1,c)
        macz(a,b-1,c)=max(v2)
      endif
    endelse

  endfor
 
  flux=fltarr(n,n,n) ;....array of log10(Flux[erg/s])
  flux(*,*,*)=-33
  macx=sqrt(macx^2.+macy^2.+macz^2.)
  i3=where(macx gt mthr,nn)

  print,nn,'shocks M>',mthr
  bo=array_indices(flux,i3)
  for i=0L,nn-1 do begin
    a=bo(0,i)
    b=bo(1,i)
    c=bo(2,i)
    m=macx(a,b,c)
    etha=0.
    delta=0.

    ; goto,ryu^M
    angamma=1.66666
    ga=1.66666
    ; rr=(ga+1.)/float(ga-1.+2./float(m*m))^M
    ; d0=(2./float(ga*(ga-1)*m*m*rr))*((2*ga*m*m-ga+1.)/float(ga+1.)-rr^ga)^M

    ; etha=0.^M
    ;m4=1/float(m*m*m*m)^M
    ;delta=d0                        ;*0.92^M
    ;delta=-4.25/float(m4)+6.42*(m-1)/(float(m4))-1.34*(m-1)*(m-1)/float(m4)+1.26*\
    ;(m-1)*(m-1)*(m-1)/float(m4)+0.275*(m-1)*(m-1)*(m-1)*(m-1)/float(m4)^M
    ; etha=5.46*(m4)-9.78*(m-1)/(float(m4))+4.17*(m-1)*(m-1)/float(m4)-0.334*(m-1)*\
    ;(m-1)*(m-1)/float(m4)+0.57*(m-1)*(m-1)*(m-1)*(m-1)/float(m4)^M

    ;..the following would be the one to use
    ; etha=5.46*(m4)-9.78*(m-1)*(m4)+4.17*(m-1)^2.*(m4)-0.334*(m-1)^3.*(m4)+0.57*(m-1)^4.*(m4)^M


    dem=dens(a,b,c)*(m*m+3)/float(4*m*m)
    tem=temp(a,b,c)*(16*m*m)/float((5*m*m-1)*(m*m+3))
    ;cr(a,b,c)=alog10(0.5*etha*dem)+0.666*vcell+3*alog10(m)+1.5*alog10(angamma*kboltz*tem*mw_pm)
    flux(a,b,c)=alog10(0.5*dem)+0.666*vcell+3*alog10(m)+1.5*alog10(angamma*kboltz*tem*mw_pm)
  endfor

end

;......NORMALIZZAZIONE DEL CAMPIONE COMPLETO CON LE SCALE DI FILTRAGGIO FISSO......

pro do_norm

print,'Inizio della normalizzazione dei profili a scale mobili'

fold_out='/home/STUDENTI/matteo.angelinelli/Output/Prof/'    ;...output folder
   cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]   ;...list of available clusters (considering z=0 snapshots)
   nc=size(cl_list)
   sc=['3','5','10','20','30']

CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine

nbr=30     ;…numero di bin per campionare il raggio (va sperimentato) 
dr=1/float(nbr)   ;,bin radiale
ro=dr*(indgen(nbr)) ;…genera coordinata radiale    ;...FRANCO: aggiunt 0.5
pratio_interp=fltarr(nc(1),nbr)    ;…array di profilo radiale per tutti i cluster
pratio_plot=fltarr(nc(1),280)
rr_plot=fltarr(280,nc(1))
rr=fltarr(nc(1),nbr)
rv_vect=fltarr(nc(1))


for s=0,4 do begin 
set_plot,'x'
window,s
plot,indgen(2),indgen(2),/nodata,xrange=[1e-3,3],yrange=[1.e-4,1.],/ylog

  for c=0,nc(1)-1 do begin
  readcol,fold_out+cl_list(c)+'_pressure_'+sc(s)+'_profile.dat',r,pturb,pth,ratio,pd,pdm,FORMAT='(F,F,F,F,F,F)',SKIPLINE=1,numline=280

  pratio_plot(c,*)=ratio
  ;rr_plot(*,c)=r


  rv2=where(r gt 0.99 and r lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(r gt 0.98 and r lt 1.02,nr) 
  rv=rv2(0)                            
endif

  rv_vect(c)=rv

  rr(c,*)=interpolate(r(0:2*rv),2*rv*ro) 
  pratio_interp(c,*)=interpolate(ratio(0:2*rv),2*rv*indgen(nbr)/float(nbr),/cubic) 


openw,1,fold_out+'/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(c)+'_pressure_30_'+sc(s)+'.dat'
for i=0,nbr-1 do begin 
printf,1,rr(c,i), pratio_interp(c,i),format='(e,e)'
endfor
close,1

oplot,rr(c,*),pratio_interp(c,*)

endfor

endfor

close,1

print,'Fine della normalizzazione dei profili'

end

;.....NORMALIZZAZIONE DEL CATALOGO FILTRATTO A R200/3.....

pro do_norm_r200_3

print,'......inizio normalizzazione a R200/3.....'

fold_out='/home/STUDENTI/matteo.angelinelli/Output/Prof/'    ;...output folder
   cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]   ;...list of available clusters (considering z=0 snapshots)
   nc=size(cl_list)

CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine

nbr=30     ;…numero di bin per campionare il raggio (va sperimentato) 
dr=1/float(nbr)   ;,bin radiale
ro=dr*(indgen(nbr)) ;…genera coordinata radiale    ;...FRANCO: aggiunt 0.5
pratio_interp=fltarr(nc(1),nbr)    ;…array di profilo radiale per tutti i cluster
pratio_plot=fltarr(nc(1),280)
rr_plot=fltarr(280,nc(1))
rr=fltarr(nc(1),nbr)
rv_vect=fltarr(nc(1))
close,1
close,2

set_plot,'x'
window,0
plot,indgen(2),indgen(2),/nodata,yrange=[1.e-4,1.],/ylog

  for c=0,nc(1)-1 do begin
  readcol,fold_out+cl_list(c)+'_scalamobile_pressure_profile.dat',r,pturb,pth,ratio,pd,pdm,FORMAT='(F,F,F,F,F,F)',SKIPLINE=1,numline=280

  pratio_plot(c,*)=ratio
  ;rr_plot(*,c)=r

  rv2=where(r gt 0.99 and r lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(r gt 0.98 and r lt 1.02,nr) 
  rv=rv2(0)                            
endif
  
  print,r
  print,rv
 
  rv_vect(c)=rv

  rr(c,*)=interpolate(r(0:2*rv),2*rv*ro) 
  pratio_interp(c,*)=interpolate(ratio(0:2*rv),2*rv*indgen(nbr)/float(nbr),/cubic) 


	openw,1,fold_out+'/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(c)+'_pressure_scalamobile.dat'
	for i=0,nbr-1 do begin 
	printf,1,rr(c,i), pratio_interp(c,i),format='(e,e)'
	endfor
	close,1

	endfor

	  pratio_mean=fltarr(5,nbr)  
	for i=0,nbr-1 do begin
	 igood=where(finite(pratio_interp(*,i)) eq 1 and pratio_interp(*,i) ne 0 )
	 pratio_mean(0:4,i)=percentiles(pratio_interp(igood,i),value=[0.1,0.25,0.5,0.75,0.9])  
	endfor

	r_perc=fltarr(nbr)
	for i=0,nbr-1 do begin  
	  r_perc(i)=median(rr(*,i))
	endfor

	openw,2,fold_out+'/Prof_inter/'+'Profile_interpolate_median_pressure_scalamobile.dat'
	for i=0,nbr-1 do begin 
	printf,2,r_perc(i), pratio_mean(2,i),format='(e,e)'
	endfor
	close,2

end

;......FIT DEL CAMPIONE COMPLETO ALLE SCALE DI FILTRAGGIO FISSE CON IL MODELLO DI NELSON E NOSTRO.....

pro nelson

print,'Inizio fit dei profili con il nostro modello e il modello di Nelson'

fold_out='/home/STUDENTI/matteo.angelinelli/Output/Prof/'    ;...output folder
fold2='/home/STUDENTI/matteo.angelinelli/Output/'
   cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]   ;...list of available clusters (considering z=0 snapshots)
   nc=size(cl_list)
   sc=['3','5','10','20','30']

CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine


nbr=30   ;...è il numero dei bin scelti per campionare i raggi. entra nel nome del file generato da un'altra procedure. valori: 10,25,30,50,100
nbr_name='30'

w=fltarr(nbr)

close,1
close,3
;openw,3,fold_out+'riassunto_Nelson_pressure_2.dat'

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'Nelson_median_fit_2.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='Alpha',/xstyle,/ystyle,/ylog
color=indgen(5)*55.
line=indgen(5)*0.
colorcluster=indgen(nc(1))*20.

for s=0,4 do begin
;s=4
rr_plot=fltarr(nbr,nc(1))
rratio_plot=fltarr(nbr,nc(1))

r=fltarr(nbr)
ratio=fltarr(nbr)

for c=0,nc(1)-1 do begin 
  readcol,fold_out+'/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(c)+'_pressure_30_'+sc(s)+'.dat',rr,rratio,format='(f,f)',/silent ;...file with single cluster

rr_plot(*,c)=rr
rratio_plot(*,c)=rratio

endfor

for i=0,nbr-1 do begin 
r(i)=median(rr_plot(i,*))
ratio(i)=median(rratio_plot(i,*))
endfor

dev=fltarr(nbr)

for i=0,nbr-1 do begin 
dev(i)=stddev(rratio_plot(i,*))
endfor

f=finite(dev)
ff=where(f ne 0)
for i=0,nbr-1 do begin 
if f(i) eq 0 then begin 
dev(i)=mean(dev(ff))
endif
endfor
  w(*)=1./dev(*)^2.

if sc(s) eq 3 then begin 
  a0=[0.5,0.5,2.9]
  a1=[0.73,0.43,0.003]
endif
if sc(s) eq 5 then begin 
  a0=[0.5,0.5,3.0]
  a1=[0.73,0.43,0.003]
endif
if sc(s) eq 10 then begin 
  a0=[0.5,0.5,1.7]
  a1=[0.7,0.25,0.02]
endif
if sc(s) eq 20 or sc(s) eq 30 then begin 
  a0=[0.5,0.5,1.7]
  a1=[0.7,0.25,0.01]
endif
if sc(s) eq 30 then begin 
  a0=[0.5,0.5,1.7]
  a1=[0.7,0.25,0.01]
endif

yfit=curvefit(r(1:nbr-1),ratio(1:nbr-1),w(1:nbr-1),a0,sigma0,function_name='gfunct',itmax=1000,chisq=chisq0,yerror=yerr0,iter=niter0,status=status0)
yfit2=curvefit(r(1:nbr-1),ratio(1:nbr-1),w(1:nbr-1),a1,sigma1,function_name='funct_nelson',itmax=1000,chisq=chisq1,yerror=yerr1,iter=niter1,status=status1)


;  printf,3,'Smoothing scale: ',sc(s)*20,format='(a,f)

;  printf,3,'Our: ',' a0 ',a0(0),'+-',sigma0(0),' a1 ',a0(1),'+-',sigma0(1),' a2 ',a0(2),'+-',sigma0(2),format='(a,a,f,a,f,a,f,a,f,a,f,a,f)
;  printf,3,'Nelson: ',' A ',a1(0),'+-',sigma1(0),' B ',a1(1),'+-',sigma1(1),' Gamma ',a1(2),'+-',sigma1(2),format='(a,a,f,a,f,a,f,a,f,a,f,a,f)

print,'Smoothing scale: ',sc(s)*20,' Convergence: ',status0,status1
print,'Our: ',' a0 ',a0(0),'+-',sigma0(0),' a1 ',a0(1),'+-',sigma0(1),' a2 ',a0(2),'+-',sigma0(2)
print,'Nelson: ',' A ',a1(0),'+-',sigma1(0),' B ',a1(1),'+-',sigma1(1),' Gamma ',a1(2),'+-',sigma1(2)


chi0=0.
chi1=0.

for i=1,nbr-1 do begin 
chi0+=(ratio(i)-yfit(i-1))^2./dev(i)^2.
chi1+=(ratio(i)-yfit2(i-1))^2./dev(i)^2.
endfor

;printf,3,'Our model chi-square: ',chi0,format='(a,f)
;printf,3,'Nelson model chi-square: ',chi1,format='(a,f)

;printf,3,'------------------------',format='(a)'


print,'Our model chi-square: ',chisq0,chi0
print,'Nelson model chi-square: ',chisq1,chi1
print,'------------------------'

readcol,fold2+'eckert_point.dat',id,m200,alfa200,alfa500,format='(A,F,F,F)',skipline=2,/silent

r200=fltarr(13)
r500=fltarr(13)
r200(*)=1.
r500(*)=0.71

;if sc(s) eq 3 then begin
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_Nelson_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='E!dturb!n/E!dtot!n',/xstyle,/ystyle,/ylog
;endif
;if sc(s) ne 3 then begin
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_Nelson_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-3,1],xtitle='r/r!d200!n',ytitle='E!dturb!n/E!dtot!n',/xstyle,/ystyle,/ylog
;endif

;n=1-0.452*(1+exp(-(r/0.841)^1.628))

;set_plot,'x'
;window,s+1
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_profile_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='Alpha',/xstyle,/ystyle,/ylog



oplot,r,ratio(1:nbr-1),linestyle=2,thick=3,col=color(s)
oplot,r,yfit,thick=3,col=color(s)
;for c=0,nc(1)-1 do begin 
;oplot,r,rratio_plot(*,c),linestyle=3,col=colorcluster(c)
;endfor
;if s eq 0 then begin 
;legend,[cl_list],linestyle=3,col=[colorcluster],/left,/top
;legend,['Data','Median','Fit'],linestyle=[3,2,0],col=0,/right,/top
;legend,['60 kpc'],col=0,/left,/bottom
;endif
;if s eq 1 then begin 
;legend,[cl_list],linestyle=3,col=[colorcluster],/right,/bottom
;legend,['Data','Median','Fit'],linestyle=[3,2,0],col=0,/left,/bottom
;legend,['100 kpc'],col=0,/left,/bottom
;endif
;if s eq 2 or s eq 3 or s eq 4 then begin 
;legend,[cl_list],linestyle=3,col=[colorcluster],/right,/bottom
;legend,['Data','Median','Fit'],linestyle=[3,2,0],col=0,/left,/bottom
;endif
;if s eq 2 then begin 
;legend,['200 kpc'],col=0,/left,/bottom
;endif
;if s eq 3 then begin 
;legend,['400 kpc'],col=0,/left,/bottom
;endif
;if s eq 4 then begin 
;legend,['600 kpc'],col=0,/left,/bottom
;endif

;device,/close


oplot,r,yfit2,thick=2,linestyle=3,col=color(s)

;oplot,r200,alfa200/100.,psym=5
;oplot,r500,alfa500/100.,psym=4
;legend,['Median','Our Fit','Nelson Fit','Nelson et al. +14'],linestyle=[2,0,3,1],position=[0.0,1.1e-3];,/center_legend,/bottom_legend;,col=[0,200,100,150],/left_legend,/top_legend
;legend,['Median','Our Fit','Nelson Fit'],linestyle=[2,0,3],/left,/bottom
;device,/close

;if sc(s) eq 3 then begin
;openw,1,fold_out+'fit_parameters_Nelson.dat'    
;printf,1,'Smooth Scale ',' a0 ',' a1 ',' a2 ',' A ',' B ',' Gamma ',format='(a,a,a,a,a,a,a)'
;printf,1,sc(s),a0(0),a0(1),a0(2),a1(0),a1(1),a1(2),format='(i,e,e,e,e,e,e)'
;endif
;if sc(s) ne 3 then begin 
;printf,1,sc(s),a0(0),a0(1),a0(2),a1(0),a1(1),a1(2),format='(i,e,e,e,e,e,e)
;endif

endfor
nhform=1-0.452*(1+exp(-(r/0.841)^1.628))
oplot,r,nhform,thick=3,linestyle=1
close,1
close,3
;legend,['Eckert et al. +18 -- R500','Eckert et al. +18 -- R200'],psym=[4,5],position=[0.0,3.e-4]
legend,['Median','Our Fit','Nelson Fit','Nelson et al. +14'],linestyle=[2,0,3,1],/left,/bottom
;legend,['Median','Fit'],linestyle=[2,0],col=0,/left_legend,/bottom_legend
legend,['60kpc','100kpc','200kpc','400kpc','600kpc'],linestyle=line,col=color,/right_legend,/bottom_legend
device,/close
end


;.....FUNZIONE DEL NOSTRO MODELLO......

pro gfunct, x,a,f,pder      ; Function used to make the fit

    f=a[0]*x^a[1]+a[2]     
    pder=[[x^a[1]],[a[0]*a[1]*x^(a[1]-1)],[x*0.+1]]
end

;......FUNZIONE DEL MODELLO DI NELSON.....

pro funct_nelson, x,a,f,pder      ; Function used in Nelson et al. +14

    f=1-a[0]*(1+exp(-(x/a[1])^a[2]))
    pder=[[-exp(-(x/a[1])^a[2])-1],[-(a[0]*a[1]*exp(-(x/a[1])^a[2])*(x/a[1])^a[2])/a[1]],[a[0]*exp(-(x/a[1])^a[2])*(x/a[1])^a[2]*alog(x/a[1])]]
end

;.....DIVISIONE IN MASSA DEL CATALOGO A TUTTE LE SCALE DI FILTRAGGIO....

pro mass

fold_out='/home/STUDENTI/matteo.angelinelli/Output/'

sc=['3','5','10','20','30']

readcol,fold_out+'catalog_high_mass.dat',IDhigh,Mhh,Mth,format='(A,F,F)',/silent
readcol,fold_out+'catalog_low_mass.dat',IDlow,Mhl,Mtl,format='(A,F,F)',/silent

nc1=size(IDhigh)
nc2=size(IDlow)
nbr=30
nc=nc1(1)+nc2(1)

r_high=fltarr(nbr,nc1(1))
r_low=fltarr(nbr,nc1(1))
ratio_high=fltarr(nbr,nc1(1))
ratio_low=fltarr(nbr,nc1(1))
w_h=fltarr(nbr)
w_l=fltarr(nbr)
ws_h=fltarr(nbr)
ws_l=fltarr(nbr)
r_scaleh=fltarr(nbr,nc1(1))
ratio_scaleh=fltarr(nbr,nc1(1))
r_scalel=fltarr(nbr,nc2(1))
ratio_scalel=fltarr(nbr,nc2(1))
color=indgen(6)*55.
line=indgen(6)*0.
close,1
close,2

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'All_ourmodel_median_2.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='Alpha',/xstyle,/ystyle,/ylog


;openw,1,fold_out+'fit_parameters_high_mass_2.dat'    
;printf,1,'Smooth Scale ',' a0 ',' err_a0 ',' a1 ',' err_a1 ',' a2 ',' err_a1 ','ChiSquare',format='(a,a,a,a,a,a,a)'
;openw,2,fold_out+'fit_parameters_low_mass_2.dat'    
;printf,2,'Smooth Scale ',' a0 ',' a1 ',' a2 ','ChiSquare',format='(a,a,a,a,a,a,a,a)'
;openw,3,fold_out+'fit_parameters_Nelson_high_mass_2.dat'    
;printf,3,'Smooth Scale ',' a0 ',' a1 ',' a2 ',' A ',' B ',' Gamma ',' ChiSquareOur ',' ChiSquareNelson ',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
;openw,4,fold_out+'fit_parameters_Nelson_low_mass_2.dat'    
;printf,4,'Smooth Scale ',' a0 ',' a1 ',' a2 ',' A ',' B ',' Gamma ',' ChiSquareOur ',' ChiSquareNelson ',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'



for c=0,nc/2-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDhigh(c)+'_pressure_scalamobile.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

r_scaleh(*,c)=rr1
ratio_scaleh(*,c)=rratio1

endfor

for c=nc/2,nc-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDlow(c-nc/2)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

r_scalel(*,c-nc/2)=rr2
ratio_scalel(*,c-nc/2)=rratio2
endfor

rs_h=fltarr(nbr)
rs_l=fltarr(nbr)
ratiosh=fltarr(nbr)
ratiosl=fltarr(nbr)

for i=0,nbr-1 do begin 
rs_h(i)=median(r_scaleh(i,*))
rs_l(i)=median(r_scalel(i,*))
ratiosh(i)=median(ratio_scaleh(i,*))
ratiosl(i)=median(ratio_scalel(i,*))
endfor

devs_high=fltarr(nbr)
devs_low=fltarr(nbr)

for i=0,nbr-1 do begin 
devs_high(i)=stddev(ratio_scaleh(i,*))
devs_low(i)=stddev(ratio_scalel(i,*))
endfor

f1=finite(devs_high)
f2=finite(devs_low)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
devs_high(i)=mean(devs_high(ff1))
endif
if f2(i) eq 0 then begin 
devs_low(i)=mean(devs_low(ff2))
endif
endfor

  ws_h(*)=1./devs_high(*)^2.
  ws_l(*)=1./devs_low(*)^2.


ahs=[0.5,0.5,1.7]
als=[0.6,0.5,0.4]

anhs=[0.1,0.1,0.1]
anls=[0.1,0.1,0.1]


 yfiths=curvefit(rs_h(1:nbr-1),ratiosh(1:nbr-1),ws_h(1:nbr-1),ahs,sigmahs,function_name='gfunct',itmax=1000,chisq=chisqhs,yerror=yerrhs,iter=niterhs,status=statushs)
 yfitls=curvefit(rs_l(1:nbr-1),ratiosl(1:nbr-1),ws_l(1:nbr-1),als,sigmals,function_name='gfunct',itmax=1000,chisq=chisqls,yerror=yerrls,iter=niterls,status=statusls)

 yfit2hs=curvefit(rs_h(1:nbr-1),ratiosh(1:nbr-1),ws_h(1:nbr-1),anhs,sigma1hs,function_name='funct_nelson',itmax=1000,chisq=chisqnhs,yerror=yerrnhs,iter=niternhs,status=statusnhs)
 yfit2ls=curvefit(rs_l(1:nbr-1),ratiosl(1:nbr-1),ws_l(1:nbr-1),anls,sigma1ls,function_name='funct_nelson',itmax=1000,chisq=chisqnls,yerror=yerrnls,iter=niternls,status=statusnls)


print,'Smoothing scale: R200/3',' Convergence: ',statushs,statusls,statusnhs,statusnls
print,'High: ',' a0 ',ahs(0),'+-',sigmahs(0),' a1 ',ahs(1),'+-',sigmahs(1),' a2 ',ahs(2),'+-',sigmahs(2)
print,'Low: ',' a0 ',als(0),'+-',sigmals(0),' a1 ',als(1),'+-',sigmals(1),' a2 ',als(2),'+-',sigmals(2)
print,'Nelson, High: ',' a0 ',anhs(0),'+-',sigma1hs(0),' a1 ',anhs(1),'+-',sigma1hs(1),' a2 ',anhs(2),'+-',sigma1hs(2)
print,'Nelson, Low: ',' a0 ',anls(0),'+-',sigma1ls(0),' a1 ',anls(1),'+-',sigma1ls(1),' a2 ',anls(2),'+-',sigma1ls(2)

chihs=0.
chils=0.

chinls=0.
chinhs=0.

for i=1,nbr-1 do begin 
chihs+=(ratiosh(i)-yfiths(i-1))^2./devs_high(i)^2.
chils+=(ratiosl(i)-yfitls(i-1))^2./devs_low(i)^2.
chinhs+=(ratiosh(i)-yfit2hs(i-1))^2./devs_high(i)^2.
chinls+=(ratiosl(i)-yfit2ls(i-1))^2./devs_low(i)^2.
endfor

print,'High chi-square - Our model: ',chihs
print,'Low chi-square - Our model: ',chils
print,'High chi-square - Nelson model: ',chinhs
print,'Low chi-square - Nelson model: ',chinls
print,'------------------------'

oplot,rs_h,ratiosh(1:nbr-1),linestyle=2,thick=2,col=color(0)
oplot,rs_l,ratiosl(1:nbr-1),linestyle=0,thick=2,col=color(0)

;oplot,rs_h,yfiths,thick=2,linestyle=0,col=color(0)
;oplot,rs_l,yfitls,thick=2,linestyle=0,col=color(0)

;oplot,rs_h,yfit2hs,thick=2,linestyle=3,col=color(0)
;oplot,rs_l,yfit2ls,thick=2,linestyle=3,col=color(0)

;printf,1,'R200/3 ',ahs(0),sigmahs(0),ahs(1),sigmahs(1),ahs(2),sigmahs(2),chihs,format='(a,e,e,e,e,e,e,e)'
;printf,2,'R200/3 ',als(0),sigmals(0),als(1),sigmals(1),als(2),sigmals(2),chils,format='(a,e,e,e,e,e,e,e)'
;printf,3,'R200/3 ',ahs(0),sigmahs(0),ahs(1),sigmahs(1),ahs(2),sigmahs(2),anhs(0),sigma1hs(0),anhs(1),sigma1hs(1),anhs(2),sigma1hs(2),chihs,chinhs,format='(a,e,e,e,e,e,e,e,e,e,e,e,e,e,e)'
;printf,4,'R200/3 ',als(0),sigmals(0),als(1),sigmals(1),als(2),sigmals(2),anls(0),sigma1ls(0),anls(1),sigma1ls(1),anls(2),sigma1ls(2),chils,chinls,format='(a,e,e,e,e,e,e,e,e,e,e,e,e,e,e)'


for s=0,4 do begin 

;s=4

for c=0,nc1(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDhigh(c)+'_pressure_30_'+sc(s)+'.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

r_high(*,c)=rr1
ratio_high(*,c)=rratio1
endfor

for c=0,nc1(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDlow(c)+'_pressure_30_'+sc(s)+'.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

r_low(*,c)=rr2
ratio_low(*,c)=rratio2
endfor

r_h=fltarr(nbr)
r_l=fltarr(nbr)
ratioh=fltarr(nbr)
ratiol=fltarr(nbr)

for i=0,nbr-1 do begin 
r_h(i)=median(r_high(i,*))
r_l(i)=median(r_low(i,*))
ratioh(i)=median(ratio_high(i,*))
ratiol(i)=median(ratio_low(i,*))
endfor

dev_high=fltarr(nbr)
dev_low=fltarr(nbr)


for i=0,nbr-1 do begin 
dev_high(i)=stddev(ratio_high(i,*))
dev_low(i)=stddev(ratio_low(i,*))
endfor

f1=finite(dev_high)
f2=finite(dev_low)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
dev_high(i)=mean(dev_high(ff1))
endif
if f2(i) eq 0 then begin 
dev_low(i)=mean(dev_low(ff2))
endif
endfor

  w_h(*)=1./dev_high(*)^2.
  w_l(*)=1./dev_low(*)^2.

if sc(s) eq 3 then begin 
  ah=[0.5,0.5,2.9]
al=[0.5,0.5,2.9]
  anh=[0.73,0.43,0.003]
  anl=[0.73,0.43,0.003]
endif
if sc(s) eq 5 then begin 
  ah=[0.5,0.5,1.9]
al=[0.5,0.5,1.9]
  anh=[0.73,0.43,0.003]
anl=[0.73,0.43,0.003]
endif
if sc(s) eq 10 then begin 
  ah=[0.5,0.5,0.6]
al=[0.5,0.5,1.3]
  anh=[0.73,0.43,0.003]
anl=[0.73,0.43,0.003]
endif
if sc(s) eq 20 then begin 
  ah=[0.5,0.5,3.21]
al=[0.5,0.7,3.3]
   anh=[0.73,0.43,0.003]
anl=[0.73,0.43,0.003]
endif
if sc(s) eq 30 then begin 
  ah=[0.5,0.5,3.1]
al=[0.5,0.5,1.7]
  anh=[0.73,0.43,0.003]
anl=[0.73,0.43,0.003]
endif


 yfith=curvefit(r_h(1:nbr-1),ratioh(1:nbr-1),w_h(1:nbr-1),ah,sigmah,function_name='gfunct',itmax=1000,chisq=chisqh,yerror=yerrh,iter=niterh,status=statush)
 yfitl=curvefit(r_l(1:nbr-1),ratiol(1:nbr-1),w_l(1:nbr-1),al,sigmal,function_name='gfunct',itmax=1000,chisq=chisql,yerror=yerrl,iter=niterl,status=statusl)

yfit2h=curvefit(r_h(1:nbr-1),ratioh(1:nbr-1),w_h(1:nbr-1),anh,sigma1h,function_name='funct_nelson',itmax=1000,chisq=chisq1h,yerror=yerr1h,iter=niter1h,status=status1h)
yfit2l=curvefit(r_l(1:nbr-1),ratiol(1:nbr-1),w_l(1:nbr-1),anl,sigma1l,function_name='funct_nelson',itmax=1000,chisq=chisq1l,yerror=yerr1l,iter=niter1l,status=status1l)


print,'Smoothing scale: ',sc(s)*20,' Convergence: ',statush,statusl,status1h,status1l
print,'High: ',' a0 ',ah(0),'+-',sigmah(0),' a1 ',ah(1),'+-',sigmah(1),' a2 ',ah(2),'+-',sigmah(2)
print,'Low: ',' a0 ',al(0),'+-',sigmal(0),' a1 ',al(1),'+-',sigmal(1),' a2 ',al(2),'+-',sigmal(2)
print,'Nelson, high: ',' A ',anh(0),'+-',sigma1h(0),' B ',anh(1),'+-',sigma1h(1),' Gamma ',anh(2),'+-',sigma1h(2)
print,'Nelson, low: ',' A ',anl(0),'+-',sigma1l(0),' B ',anl(1),'+-',sigma1l(1),' Gamma ',anl(2),'+-',sigma1l(2)


chih=0.
chil=0.

chinh=0.
chinl=0.

for i=1,nbr-1 do begin 
chih+=(ratioh(i)-yfith(i-1))^2./dev_high(i)^2.
chil+=(ratiol(i)-yfitl(i-1))^2./dev_low(i)^2.
chinh+=(ratioh(i)-yfit2h(i-1))^2./dev_high(i)^2.
chinl+=(ratiol(i)-yfit2l(i-1))^2./dev_low(i)^2.
endfor

print,'High chi-square: ',chih
print,'Low chi-square: ',chil
print,'High chi-square - Nelson model: ',chinh
print,'Low chi-square - Nelson model: ',chinl
print,'------------------------'


readcol,fold_out+'eckert_point.dat',id,m200,alfa200,alfa500,format='(A,F,F,F)',skipline=2,/silent
r200=fltarr(13)
r500=fltarr(13)
r200(*)=1.
r500(*)=0.71

med_m=median(m200)
sm=size(m200)

a200_low=fltarr(sm(1)/2+1)
a200_high=fltarr(sm(1)/2)

a500_low=fltarr(sm(1)/2+1)
a500_high=fltarr(sm(1)/2)

i=0
j=0

for c=0,sm(1)-1 do begin 
if m200(c) le med_m then begin 
a200_low(i)=alfa200(c)
a500_low(i)=alfa500(c)
i=i+1
endif
if m200(c) gt med_m then begin 
a200_high(j)=alfa200(c)
a500_high(j)=alfa500(c)
j=j+1
endif
endfor

;set_plot,'x'
;window,s

;if sc(s) eq 3 then begin
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_Nelson_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='E!dturb!n/E!dtot!n',/xstyle,/ystyle,/ylog
;endif
;if sc(s) ne 3 then begin
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_Nelson_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-3,1],xtitle='r/r!d200!n',ytitle='E!dturb!n/E!dtot!n',/xstyle,/ystyle,/ylog
;endif

;nhform=1-0.452*(1+exp(-(r_h/0.841)^1.628))
;nlform=1-0.452*(1+exp(-(r_l/0.841)^1.628))

oplot,r_h,ratioh(1:nbr-1),linestyle=2,thick=2,col=color(s+1)
oplot,r_l,ratiol(1:nbr-1),linestyle=0,thick=2,col=color(s+1)

;oplot,r_h,yfith,thick=2,linestyle=0,col=color(s+1)
;oplot,r_l,yfitl,thick=2,linestyle=0,col=color(s+1)

;oplot,r_h,yfit2h,thick=2,linestyle=3,col=color(s+1)
;oplot,r_l,yfit2l,thick=2,linestyle=3,col=color(s+1)

;oplot,r_h,nhform,linestyle=1,thick=3
;oplot,r_l,nlform,linestyle=1,thick=3

;oplot,r200,a200_high/100.,psym=5
;oplot,r500,a500_high/100.,psym=4
;oplot,r200,a200_low/100.,psym=5
;oplot,r500,a500_low/100.,psym=4



;printf,1,sc(s),ah(0),sigmah(0),ah(1),sigmah(1),ah(2),sigmah(2),chih,format='(i,e,e,e,e,e,e,e)'
;printf,2,sc(s),al(0),sigmal(0),al(1),sigmal(1),al(2),sigmal(2),chil,format='(i,e,e,e,e,e,e,e)'
;printf,3,sc(s),ah(0),sigmah(0),ah(1),sigmah(1),ah(2),sigmah(2),anh(0),sigma1h(0),anh(1),sigma1h(1),anh(2),sigma1h(2),chih,chinh,format='(i,e,e,e,e,e,e,e,e,e,e,e,e,e,e)'
;printf,4,sc(s),al(0),sigmal(0),al(1),sigmal(1),al(2),sigmal(2),anl(0),sigma1l(0),anl(1),sigma1l(1),anl(2),sigma1l(2),chil,chinl,format='(i,e,e,e,e,e,e,e,e,e,e,e,e,e,e)'



endfor

close,1
close,2
close,3
close,4
;legend,['Median','Our Fit','Nelson Fit','Nelson et al. +14'],linestyle=[2,0,3,1],/left,/bottom
;legend,['Median','Our Fit'],linestyle=[2,0],/left,/bottom
;legend,['Eckert et al. +18 -- R500','Eckert et al. +18 -- R200'],psym=[4,5],position=[1.11,2.2e-3]
;legend,['Median Low','Median High','Fit Low','Fit High'],linestyle=[2,4,0,3],/left,/bottom
legend,['Less Massive','More Massive'],linestyle=[0,2],/left,/bottom
legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],linestyle=line,col=color,/right_legend,/bottom_legend
device,/close

end

;.....ANALISI DEL CAMPIONE A TUTTE LE SCALE CON IL PARAMETRO C

pro perturbc

fold_out='/home/STUDENTI/matteo.angelinelli/Output/'
fold_file='/home/STUDENTI/matteo.angelinelli/Procedure/'

sc=['3','5','10','20','30']

readcol,fold_file+'param_sort_mean.dat',IDc,IDw,format='(A,A)',/silent,numline=14,skipline=2

nc=size(IDc)
nw=size(IDw)
nbr=30

rc_high=fltarr(nbr,nc(1)/2)
rc_low=fltarr(nbr,nc(1)/2)
rc_scaleh=fltarr(nbr,nc(1)/2)
rc_scalel=fltarr(nbr,nc(1)/2)

ratioc_high=fltarr(nbr,nc(1)/2)
ratioc_low=fltarr(nbr,nc(1)/2)
ratioc_scaleh=fltarr(nbr,nc(1)/2)
ratioc_scalel=fltarr(nbr,nc(1)/2)

wc_h=fltarr(nbr)
wc_l=fltarr(nbr)
ws_h=fltarr(nbr)
ws_l=fltarr(nbr)

color=indgen(6)*55.
line=indgen(6)*0.

close,1
close,2

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'Zoom_c_log.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,0.5],yrange=[4e-4,2e-1],xtitle='r/r!d200!n',ytitle='Alpha',/xstyle,/ystyle,/ylog

;openw,1,fold_out+'fit_parameters_high_perturb_c_2.dat'    
;printf,1,'Smooth Scale ',' a0 ',' err_a0 ',' a1 ',' err_a1 ',' a2 ',' err_a1 ','ChiSquare',format='(a,a,a,a,a,a,a,a)'
;openw,2,fold_out+'fit_parameters_low_perturb_c_2.dat'    
;printf,2,'Smooth Scale ',' a0 ',' err_a0 ',' a1 ',' err_a1 ',' a2 ',' err_a1 ','ChiSquare',format='(a,a,a,a,a,a,a,a)'

for c=0,nc(1)/2-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDc(c)+'_pressure_scalamobile.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

rc_scaleh(*,c)=rr1
ratioc_scaleh(*,c)=rratio1
endfor

for c=nc(1)/2,nc(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDc(c)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

rc_scalel(*,c-nc(1)/2)=rr2
ratioc_scalel(*,c-nc(1)/2)=rratio2
endfor

rs_h=fltarr(nbr)
rs_l=fltarr(nbr)
ratiosh=fltarr(nbr)
ratiosl=fltarr(nbr)

for i=0,nbr-1 do begin 
rs_h(i)=median(rc_scaleh(i,*))
rs_l(i)=median(rc_scalel(i,*))
ratiosh(i)=median(ratioc_scaleh(i,*))
ratiosl(i)=median(ratioc_scalel(i,*))
endfor

devs_high=fltarr(nbr)
devs_low=fltarr(nbr)

for i=0,nbr-1 do begin 
devs_high(i)=stddev(ratioc_scaleh(i,*))
devs_low(i)=stddev(ratioc_scalel(i,*))
endfor

f1=finite(devs_high)
f2=finite(devs_low)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
devs_high(i)=mean(devs_high(ff1))
endif
if f2(i) eq 0 then begin 
devs_low(i)=mean(devs_low(ff2))
endif
endfor

  ws_h(*)=1./devs_high(*)^2.
  ws_l(*)=1./devs_low(*)^2.


ahs=[0.6,0.5,0.5]
als=[0.5,0.5,1.5]


 yfiths=curvefit(rs_h(1:nbr-1),ratiosh(1:nbr-1),ws_h(1:nbr-1),ahs,sigmahs,function_name='gfunct',itmax=1000,chisq=chisqhs,yerror=yerrhs,iter=niterhs,status=statushs)
 yfitls=curvefit(rs_l(1:nbr-1),ratiosl(1:nbr-1),ws_l(1:nbr-1),als,sigmals,function_name='gfunct',itmax=1000,chisq=chisqls,yerror=yerrls,iter=niterls,status=statusls)


print,'Smoothing scale: R200/3',' Convergence: ',statushs,statusls
print,'High: ',' a0 ',ahs(0),'+-',sigmahs(0),' a1 ',ahs(1),'+-',sigmahs(1),' a2 ',ahs(2),'+-',sigmahs(2)
print,'Low: ',' a0 ',als(0),'+-',sigmals(0),' a1 ',als(1),'+-',sigmals(1),' a2 ',als(2),'+-',sigmals(2)

chihs=0.
chils=0.

for i=1,nbr-1 do begin 
chihs+=(ratiosh(i)-yfiths(i-1))^2./devs_high(i)^2.
chils+=(ratiosl(i)-yfitls(i-1))^2./devs_low(i)^2.
endfor

print,'High chi-square: ',chihs
print,'Low chi-square: ',chils
print,'------------------------'

oplot,rs_h,ratiosh(1:nbr-1),linestyle=2,thick=2,col=color(0)
oplot,rs_l,ratiosl(1:nbr-1),linestyle=0,thick=2,col=color(0)

;oplot,rs_h,yfiths,thick=2,linestyle=0,col=color(0)
;oplot,rs_l,yfitls,thick=2,linestyle=0,col=color(0)

;printf,1,'R200/3 ',ahs(0),sigmahs(0),ahs(1),sigmahs(1),ahs(2),sigmahs(2),chihs,format='(a,e,e,e,e,e,e,e)'
;printf,2,'R200/3 ',als(0),sigmals(0),als(1),sigmals(1),als(2),sigmals(2),chils,format='(a,e,e,e,e,e,e,e)'


for s=0,4 do begin 
;s=4

for c=0,nc(1)/2-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDc(c)+'_pressure_30_'+sc(s)+'.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

rc_high(*,c)=rr1
ratioc_high(*,c)=rratio1
endfor

for c=nc(1)/2,nc(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDc(c)+'_pressure_30_'+sc(s)+'.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

rc_low(*,c-nc(1)/2)=rr2
ratioc_low(*,c-nc(1)/2)=rratio2
endfor

rc_h=fltarr(nbr)
rc_l=fltarr(nbr)
ratioch=fltarr(nbr)
ratiocl=fltarr(nbr)

for i=0,nbr-1 do begin 
rc_h(i)=median(rc_high(i,*))
rc_l(i)=median(rc_low(i,*))
ratioch(i)=median(ratioc_high(i,*))
ratiocl(i)=median(ratioc_low(i,*))
endfor

devc_high=fltarr(nbr)
devc_low=fltarr(nbr)


for i=0,nbr-1 do begin 
devc_high(i)=stddev(ratioc_high(i,*))
devc_low(i)=stddev(ratioc_low(i,*))
endfor

f1=finite(devc_high)
f2=finite(devc_low)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
devc_high(i)=mean(devc_high(ff1))
endif
if f2(i) eq 0 then begin 
devc_low(i)=mean(devc_low(ff2))
endif
endfor

  wc_h(*)=1./devc_high(*)^2.
  wc_l(*)=1./devc_low(*)^2.

if sc(s) eq 3 then begin 
  ahc=[0.6,0.4,0.1]
alc=[0.6,0.4,0.9]

endif
if sc(s) eq 5 then begin 
  ahc=[0.5,0.5,2.1]
alc=[0.5,0.5,1.1]
  
endif
if sc(s) eq 10 then begin 
  ahc=[0.6,0.5,1.7]
alc=[0.5,0.5,0.7]

endif
if sc(s) eq 20 then begin 
  ahc=[0.5,0.5,2.35]
alc=[0.6,0.4,0.6]

endif
if sc(s) eq 30 then begin 
  ahc=[0.6,0.4,0.8]
alc=[0.5,0.5,1.5]

endif


 yfithc=curvefit(rc_h(1:nbr-1),ratioch(1:nbr-1),wc_h(1:nbr-1),ahc,sigmach,function_name='gfunct',itmax=1000,chisq=chisqhc,yerror=yerrhc,iter=niterhc,status=statushc)
 yfitlc=curvefit(rc_l(1:nbr-1),ratiocl(1:nbr-1),wc_l(1:nbr-1),alc,sigmacl,function_name='gfunct',itmax=1000,chisq=chisqlc,yerror=yerrlc,iter=niterlc,status=statuslc)

print,'Smoothing scale: ',sc(s)*20,' Convergence: ',statushc,statuslc
print,'High: ',' a0 ',ahc(0),'+-',sigmach(0),' a1 ',ahc(1),'+-',sigmach(1),' a2 ',ahc(2),'+-',sigmach(2)
print,'Low: ',' a0 ',alc(0),'+-',sigmacl(0),' a1 ',alc(1),'+-',sigmacl(1),' a2 ',alc(2),'+-',sigmacl(2)

chihc=0.
chilc=0.



for i=1,nbr-1 do begin 
chihc+=(ratioch(i)-yfithc(i-1))^2./devc_high(i)^2.
chilc+=(ratiocl(i)-yfitlc(i-1))^2./devc_low(i)^2.

endfor

print,'High chi-square: ',chihc
print,'Low chi-square: ',chilc

print,'------------------------'


readcol,fold_file+'eckert_red_error_K0.dat',id,m200,alfa200,alfa200p,alfa200m,alfa500,alfa500p,alfa500m,k0,format='(A,F,F,F,F,F,F,F,F)',skipline=2,/silent
r200=fltarr(13)
r500=fltarr(13)
r200(*)=1.
r500(*)=0.71

med_k=median(k0)
sk=size(k0)

a200_low=fltarr(sk(1)/2+1)
a200_high=fltarr(sk(1)/2)

a500_low=fltarr(sk(1)/2+1)
a500_high=fltarr(sk(1)/2)

i=0
j=0

for c=0,sk(1)-1 do begin 
if k0(c) le med_k then begin 
a200_low(i)=alfa200(c)
a500_low(i)=alfa500(c)
i=i+1
endif
if k0(c) gt med_k then begin 
a200_high(j)=alfa200(c)
a500_high(j)=alfa500(c)
j=j+1
endif
endfor


oplot,rc_h,ratioch(1:nbr-1),linestyle=2,thick=2,col=color(s+1)
oplot,rc_l,ratiocl(1:nbr-1),linestyle=0,thick=2,col=color(s+1)

;oplot,rc_h,yfithc,thick=2,linestyle=0,col=color(s+1)
;oplot,rc_l,yfitlc,thick=2,linestyle=0,col=color(s+1)

;oplot,r200,a200_high/100.,psym=5
;oplot,r500,a500_high/100.,psym=4
;oplot,r200,a200_low/100.,psym=5
;oplot,r500,a500_low/100.,psym=4

;printf,1,sc(s),ahc(0),sigmach(0),ahc(1),sigmach(1),ahc(2),sigmach(2),chihc,format='(i,e,e,e,e,e,e,e)'
;printf,2,sc(s),alc(0),sigmacl(0),alc(1),sigmacl(1),alc(2),sigmacl(2),chilc,format='(i,e,e,e,e,e,e,e)'



endfor

close,1
close,2

;legend,['Median','Our Fit'],linestyle=[2,0],/left_legend,/bottom_legend
;legend,['Eckert et al. +18 -- R500','Eckert et al. +18 -- R200'],psym=[4,5],position=[0.57,3.e-4]
legend,['Relaxed','Perturbed'],linestyle=[0,2],/left_legend,/bottom_legend
legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],linestyle=line,col=color,/right_legend,/bottom_legend
device,/close

end

;.....ANALISI DEL CAMPIONE A TUTTE LE SCALE CON IL PARAMETRO W.....

pro perturbw

fold_out='/home/STUDENTI/matteo.angelinelli/Output/'
fold_file='/home/STUDENTI/matteo.angelinelli/Procedure/'

sc=['3','5','10','20','30']

readcol,fold_file+'param_sort_mean.dat',IDc,IDw,format='(A,A)',/silent,numline=14,skipline=2

nc=size(IDc)
nw=size(IDw)
nbr=30

rw_high=fltarr(nbr,nw(1)/2)
rw_low=fltarr(nbr,nw(1)/2)
rw_scaleh=fltarr(nbr,nw(1)/2)
rw_scalel=fltarr(nbr,nw(1)/2)

ratiow_high=fltarr(nbr,nw(1)/2)
ratiow_low=fltarr(nbr,nw(1)/2)
ratiow_scaleh=fltarr(nbr,nw(1)/2)
ratiow_scalel=fltarr(nbr,nw(1)/2)

ww_h=fltarr(nbr)
ww_l=fltarr(nbr)
ws_h=fltarr(nbr)
ws_l=fltarr(nbr)

color=indgen(6)*55.
line=indgen(6)*0.

close,1
close,2

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'All_perturb_w_ourmodel_2.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='Alpha',/xstyle,/ystyle,/ylog

;openw,1,fold_out+'fit_parameters_high_perturb_w_2.dat'    
;printf,1,'Smooth Scale ',' a0 ',' err_a0 ',' a1 ',' err_a1 ',' a2 ',' err_a1 ','ChiSquare',format='(a,a,a,a,a,a,a,a)'
;openw,2,fold_out+'fit_parameters_low_perturb_w_2.dat'    
;printf,2,'Smooth Scale ',' a0 ',' err_a0 ',' a1 ',' err_a1 ',' a2 ',' err_a1 ','ChiSquare',format='(a,a,a,a,a,a,a,a)'

for c=0,nw(1)/2-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDw(c)+'_pressure_scalamobile.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

rw_scaleh(*,c)=rr1
ratiow_scaleh(*,c)=rratio1
endfor

for c=nw(1)/2,nw(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDw(c)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

rw_scalel(*,c-nw(1)/2)=rr2
ratiow_scalel(*,c-nw(1)/2)=rratio2
endfor

rs_h=fltarr(nbr)
rs_l=fltarr(nbr)
ratiosh=fltarr(nbr)
ratiosl=fltarr(nbr)

for i=0,nbr-1 do begin 
rs_h(i)=median(rw_scaleh(i,*))
rs_l(i)=median(rw_scalel(i,*))
ratiosh(i)=median(ratiow_scaleh(i,*))
ratiosl(i)=median(ratiow_scalel(i,*))
endfor

devs_high=fltarr(nbr)
devs_low=fltarr(nbr)

for i=0,nbr-1 do begin 
devs_high(i)=stddev(ratiow_scaleh(i,*))
devs_low(i)=stddev(ratiow_scalel(i,*))
endfor

f1=finite(devs_high)
f2=finite(devs_low)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
devs_high(i)=mean(devs_high(ff1))
endif
if f2(i) eq 0 then begin 
devs_low(i)=mean(devs_low(ff2))
endif
endfor

  ws_h(*)=1./devs_high(*)^2.
  ws_l(*)=1./devs_low(*)^2.


ahs=[0.6,0.4,0.6]
als=[0.6,0.4,0.8]


 yfiths=curvefit(rs_h(1:nbr-1),ratiosh(1:nbr-1),ws_h(1:nbr-1),ahs,sigmahs,function_name='gfunct',itmax=1000,chisq=chisqhs,yerror=yerrhs,iter=niterhs,status=statushs)
 yfitls=curvefit(rs_l(1:nbr-1),ratiosl(1:nbr-1),ws_l(1:nbr-1),als,sigmals,function_name='gfunct',itmax=1000,chisq=chisqls,yerror=yerrls,iter=niterls,status=statusls)


print,'Smoothing scale: R200/3',' Convergence: ',statushs,statusls
print,'High: ',' a0 ',ahs(0),'+-',sigmahs(0),' a1 ',ahs(1),'+-',sigmahs(1),' a2 ',ahs(2),'+-',sigmahs(2)
print,'Low: ',' a0 ',als(0),'+-',sigmals(0),' a1 ',als(1),'+-',sigmals(1),' a2 ',als(2),'+-',sigmals(2)

chihs=0.
chils=0.

for i=1,nbr-1 do begin 
chihs+=(ratiosh(i)-yfiths(i-1))^2./devs_high(i)^2.
chils+=(ratiosl(i)-yfitls(i-1))^2./devs_low(i)^2.
endfor

print,'High chi-square: ',chihs
print,'Low chi-square: ',chils
print,'------------------------'

oplot,rs_h,ratiosh(1:nbr-1),linestyle=2,thick=2,col=color(0)
oplot,rs_l,ratiosl(1:nbr-1),linestyle=0,thick=2,col=color(0)

;oplot,rs_h,yfiths,thick=2,linestyle=0,col=color(0)
;oplot,rs_l,yfitls,thick=2,linestyle=0,col=color(0)

;printf,1,'R200/3 ',ahs(0),sigmahs(0),ahs(1),sigmahs(1),ahs(2),sigmahs(2),chihs,format='(a,e,e,e,e,e,e,e)'
;printf,2,'R200/3 ',als(0),sigmals(0),als(1),sigmals(1),als(2),sigmals(2),chils,format='(a,e,e,e,e,e,e,e)'


for s=0,4 do begin 

;s=4

for c=0,nw(1)/2-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDw(c)+'_pressure_30_'+sc(s)+'.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

rw_high(*,c)=rr1
ratiow_high(*,c)=rratio1
endfor

for c=nw(1)/2,nw(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDw(c)+'_pressure_30_'+sc(s)+'.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

rw_low(*,c-nw(1)/2)=rr2
ratiow_low(*,c-nw(1)/2)=rratio2
endfor

rw_h=fltarr(nbr)
rw_l=fltarr(nbr)
ratiowh=fltarr(nbr)
ratiowl=fltarr(nbr)

for i=0,nbr-1 do begin 
rw_h(i)=median(rw_high(i,*))
rw_l(i)=median(rw_low(i,*))
ratiowh(i)=median(ratiow_high(i,*))
ratiowl(i)=median(ratiow_low(i,*))
endfor

devw_high=fltarr(nbr)
devw_low=fltarr(nbr)


for i=0,nbr-1 do begin 
devw_high(i)=stddev(ratiow_high(i,*))
devw_low(i)=stddev(ratiow_low(i,*))
endfor

f1=finite(devw_high)
f2=finite(devw_low)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
devw_high(i)=mean(devw_high(ff1))
endif
if f2(i) eq 0 then begin 
devw_low(i)=mean(devw_low(ff2))
endif
endfor

  ww_h(*)=1./devw_high(*)^2.
  ww_l(*)=1./devw_low(*)^2.

if sc(s) eq 3 then begin 
  ahw=[0.5,0.5,3.6]
alw=[0.5,0.5,1.6]

endif
if sc(s) eq 5 then begin 
  ahw=[0.7,0.7,1.9]
alw=[0.7,0.7,1.9]
  
endif
if sc(s) eq 10 then begin 
  ahw=[0.7,0.7,1.5]
alw=[0.6,0.6,1.5]

endif
if sc(s) eq 20 then begin 
  ahw=[0.6,0.6,2.2]
alw=[0.7,0.7,2.9]

endif
if sc(s) eq 30 then begin 
  ahw=[0.6,0.5,2.5]
alw=[0.5,0.5,3.3]

endif


 yfithw=curvefit(rw_h(1:nbr-1),ratiowh(1:nbr-1),ww_h(1:nbr-1),ahw,sigmawh,function_name='gfunct',itmax=1000,chisq=chisqhw,yerror=yerrhw,iter=niterhw,status=statushw)
 yfitlw=curvefit(rw_l(1:nbr-1),ratiowl(1:nbr-1),ww_l(1:nbr-1),alw,sigmawl,function_name='gfunct',itmax=1000,chisq=chisqlw,yerror=yerrlw,iter=niterlw,status=statuslw)

print,'Smoothing scale: ',sc(s)*20,' Convergence: ',statushw,statuslw
print,'High: ',' a0 ',ahw(0),'+-',sigmawh(0),' a1 ',ahw(1),'+-',sigmawh(1),' a2 ',ahw(2),'+-',sigmawh(2)
print,'Low: ',' a0 ',alw(0),'+-',sigmawl(0),' a1 ',alw(1),'+-',sigmawl(1),' a2 ',alw(2),'+-',sigmawl(2)

chihw=0.
chilw=0.

for i=1,nbr-1 do begin 
chihw+=(ratiowh(i)-yfithw(i-1))^2./devw_high(i)^2.
chilw+=(ratiowl(i)-yfitlw(i-1))^2./devw_low(i)^2.

endfor

print,'High chi-square: ',chihw
print,'Low chi-square: ',chilw

print,'------------------------'


readcol,fold_file+'eckert_red_error_K0.dat',id,m200,alfa200,alfa200p,alfa200m,alfa500,alfa500p,alfa500m,k0,format='(A,F,F,F,F,F,F,F,F)',skipline=2,/silent
r200=fltarr(13)
r500=fltarr(13)
r200(*)=1.
r500(*)=0.71

med_k=median(k0)
sk=size(k0)

a200_low=fltarr(sk(1)/2+1)
a200_high=fltarr(sk(1)/2)

a500_low=fltarr(sk(1)/2+1)
a500_high=fltarr(sk(1)/2)

i=0
j=0

for c=0,sk(1)-1 do begin 
if k0(c) le med_k then begin 
a200_low(i)=alfa200(c)
a500_low(i)=alfa500(c)
i=i+1
endif
if k0(c) gt med_k then begin 
a200_high(j)=alfa200(c)
a500_high(j)=alfa500(c)
j=j+1
endif
endfor

;set_plot,'x'
;window,s

;if sc(s) eq 3 then begin
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_Nelson_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1],xtitle='r/r!d200!n',ytitle='E!dturb!n/E!dtot!n',/xstyle,/ystyle,/ylog
;endif
;if sc(s) ne 3 then begin
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'Plot_Nelson_'+sc(s)+'.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-3,1],xtitle='r/r!d200!n',ytitle='E!dturb!n/E!dtot!n',/xstyle,/ystyle,/ylog
;endif

oplot,rw_h,ratiowh(1:nbr-1),linestyle=2,thick=2,col=color(s+1)
oplot,rw_l,ratiowl(1:nbr-1),linestyle=0,thick=2,col=color(s+1)

;oplot,rw_h,yfithw,thick=2,linestyle=0,col=color(s+1)
;oplot,rw_l,yfitlw,thick=2,linestyle=0,col=color(s+1)

;oplot,r200,a200_high/100.,psym=5
;oplot,r500,a500_high/100.,psym=4
;oplot,r200,a200_low/100.,psym=5
;oplot,r500,a500_low/100.,psym=4
 
;printf,1,sc(s),ahw(0),sigmawh(0),ahw(1),sigmawh(1),ahw(2),sigmawh(2),chihw,format='(i,e,e,e,e,e,e,e)'
;printf,2,sc(s),alw(0),sigmawl(0),alw(1),sigmawl(1),alw(2),sigmawl(2),chilw,format='(i,e,e,e,e,e,e,e)'

endfor

close,1
close,2

;legend,['Median','Our Fit'],linestyle=[2,0],/left_legend,/bottom_legend
;legend,['Eckert et al. +18 -- R500','Eckert et al. +18 -- R200'],psym=[4,5],position=[0.5,3.e-4]
legend,['Relaxed','Perturbed'],linestyle=[0,2],/left_legend,/bottom_legend
legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],linestyle=line,col=color,/right_legend,/bottom_legend
device,/close

end

;.....DIVISIONE IN MASSA PER IL CATALOGO FILTRATTO A R200/3........

pro mass_r200_3

fold_out='/home/STUDENTI/matteo.angelinelli/Output/'

readcol,fold_out+'catalog_high_mass.dat',IDhigh,Mhh,Mth,format='(A,F,F)',/silent
readcol,fold_out+'catalog_low_mass.dat',IDlow,Mhl,Mtl,format='(A,F,F)',/silent

nc1=size(IDhigh)
nc2=size(IDlow)
nc=nc1(1)+nc2(1)
nbr=30

r_high=fltarr(nbr,nc1(1))
r_low=fltarr(nbr,nc1(1))
ratio_high=fltarr(nbr,nc1(1))
ratio_low=fltarr(nbr,nc1(1))
w_h=fltarr(nbr)
w_l=fltarr(nbr)
w_a=fltarr(nbr)
color=indgen(3)*55.
line=indgen(2)*0.
close,1
close,2

set_plot,'x'
window,0
;set_plot,'ps'
;loadct,13
;device,filename=fold_out+'R200_high_low_mass_fit.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-2,1],xtitle='r/r!d200!n',ytitle='P!dnt!n/P!dtot!n',/xstyle,/ystyle,/ylog

readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_median_pressure_scalamobile.dat',rr_mean,rratio_mean,format='(f,f)',/silent

for c=0,nc1(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDhigh(c)+'_pressure_scalamobile.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

r_high(*,c)=rr1
ratio_high(*,c)=rratio1
endfor

for c=0,nc2(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+IDlow(c)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

r_low(*,c)=rr2
ratio_low(*,c)=rratio2
endfor

ratio_all=fltarr(nbr,nc)
for c=0,nc1(1)-1 do begin 
ratio_all(*,c)=ratio_low(*,c)
endfor
for c=0,nc2(1)-1 do begin 
ratio_all(*,c+nc1(1))=ratio_high(*,c)
endfor

r_h=fltarr(nbr)
r_l=fltarr(nbr)
ratioh=fltarr(nbr)
ratiol=fltarr(nbr)

for i=0,nbr-1 do begin 
r_h(i)=median(r_high(i,*))
r_l(i)=median(r_low(i,*))
ratioh(i)=median(ratio_high(i,*))
ratiol(i)=median(ratio_low(i,*))
endfor

dev_high=fltarr(nbr)
dev_low=fltarr(nbr)
dev_all=fltarr(nbr)

for i=0,nbr-1 do begin 
dev_high(i)=stddev(ratio_high(i,*))
dev_low(i)=stddev(ratio_low(i,*))
dev_all(i)=stddev(ratio_all(i,*))
endfor

f1=finite(dev_high)
f2=finite(dev_low)
f3=finite(dev_all)
ff1=where(f1 ne 0)
ff2=where(f2 ne 0)
ff3=where(f3 ne 0)
for i=0,nbr-1 do begin 
if f1(i) eq 0 then begin 
dev_high(i)=mean(dev_high(ff1))
endif
if f2(i) eq 0 then begin 
dev_low(i)=mean(dev_low(ff2))
endif
if f3(i) eq 0 then begin 
dev_all(i)=mean(dev_all(ff3))
endif
endfor

  w_h(*)=1./dev_high(*)^2.
  w_l(*)=1./dev_low(*)^2.
  w_a(*)=1./dev_all(*)^2.

ah=[0.1,0.1,0.1]
al=[0.1,0.2,0.5]
aa=[0.1,0.1,0.4]

 yfith=curvefit(r_h(1:nbr-1),ratioh(1:nbr-1),w_h(1:nbr-1),ah,sigmah,function_name='gfunct',itmax=1000,chisq=chisqh,yerror=yerrh,iter=niterh,status=statush)
 yfitl=curvefit(r_l(1:nbr-1),ratiol(1:nbr-1),w_l(1:nbr-1),al,sigmal,function_name='gfunct',itmax=1000,chisq=chisql,yerror=yerrl,iter=niterl,status=statusl)
 yfita=curvefit(rr_mean(1:nbr-1),rratio_mean(1:nbr-1),w_a(1:nbr-1),aa,sigmaa,function_name='gfunct',itmax=1000,chisq=chisqa,yerror=yerra,iter=nitera,status=statusa)

print,'Smoothing scale: R200/3',' Convergence: ',statush,statusl,statusa
print,'High: ',' a0 ',ah(0),'+-',sigmah(0),' a1 ',ah(1),'+-',sigmah(1),' a2 ',ah(2),'+-',sigmah(2)
print,'Low: ',' a0 ',al(0),'+-',sigmal(0),' a1 ',al(1),'+-',sigmal(1),' a2 ',al(2),'+-',sigmal(2)
print,'All: ',' a0 ',aa(0),'+-',sigmaa(0),' a1 ',aa(1),'+-',sigmaa(1),' a2 ',aa(2),'+-',sigmaa(2)

chih=0.
chil=0.
chia=0.

for i=1,nbr-1 do begin 
chih+=(ratioh(i)-yfith(i-1))^2./dev_high(i)^2.
chil+=(ratiol(i)-yfitl(i-1))^2./dev_low(i)^2.
chia+=(rratio_mean(i)-yfita(i-1))^2./dev_all(i)^2.
endfor

print,'High chi-square: ',chih
print,'Low chi-square: ',chil
print,'All chi-square: ',chia
print,'------------------------'

close,1
;openw,1,fold_out+'fit_parameters_R200.dat'    
;printf,1,'Sample',' a0 ',' err_a0 ',' a1 ',' err_a1 ',' a2 ',' err_a2 ','ChiSquare',format='(a,a,a,a,a,a,a,a)'
;printf,1,'High',ah(0),sigmah(0),ah(1),sigmah(1),ah(2),sigmah(2),chih,format='(a,e,e,e,e,e,e,e)'
;printf,1,'Low',al(0),sigmal(0),al(1),sigmal(1),al(2),sigmal(2),chil,format='(a,e,e,e,e,e,e,e)'
;printf,1,'All',aa(0),sigmaa(0),aa(1),sigmaa(1),aa(2),sigmaa(2),chia,format='(a,e,e,e,e,e,e,e)'
close,1

oplot,r_h,ratioh(1:nbr-1),linestyle=2,thick=2,col=color(1)
oplot,r_l,ratiol(1:nbr-1),linestyle=0,thick=2,col=color(2)

oplot,r_h,yfith,thick=2,linestyle=2,col=color(1)
oplot,r_l,yfitl,thick=2,linestyle=0,col=color(2)

oplot,rr_mean,rratio_mean,thick=2,linestyle=0
oplot,rr_mean,yfita,thick=2,linestyle=2

legend,['Median Low','Median High','Fit Low','Fit High'],linestyle=[0,2,0,2],col=[color(2),color(1),color(2),color(1) ],/center_legend,/bottom_legend
;legend,['Median all cluster','Fit all cluster'],linestyle=[0,2],/center_legend,/bottom_legend


;device,/close

end

;.....ANALISI DELLA RELAZIONE ALPHA-PARAMETRI X.....

pro paraxalpha

sc=['3','5','10','20','30']
fold_out='/home/STUDENTI/matteo.angelinelli/Procedure/'
fold2='/home/STUDENTI/matteo.angelinelli/Output/'
readcol,fold_out+'paraxmedi.dat',ID,cmean,wmean,format='(a,f,f)',skipline=1,/silent

n=size(ID)

ratio=fltarr(n(1),5)
ratios=fltarr(n(1))

for s=0,4 do begin 
for c=0,n(1)-1 do begin 
  readcol,fold2+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+ID(c)+'_pressure_30_'+sc(s)+'.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

  rv2=where(rr1 gt 0.99 and rr1 lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(rr1 gt 0.98 and rr1 lt 1.02,nr) 
  rv=rv2(0)                            
endif
 
  ratio(c,s)=rratio1(rv)

endfor
endfor

for c=0,n(1)-1 do begin 
  readcol,fold2+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+ID(c)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

  rv2=where(rr2 gt 0.99 and rr2 lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(rr2 gt 0.98 and rr2 lt 1.02,nr) 
  rv=rv2(0)                            
endif

  ratios(c)=rratio2(rv)

endfor

color=indgen(6)*55.
sim=[4,5,6,4,5,6]
;set_plot,'x'
;window,0
;set_plot,'ps'
;loadct,13
;device,filename=fold2+'c-alpha_2.eps',/color
;plot,indgen(2),indgen(2),/nodata,xrange=[1e-3,0.2],yrange=[1e-2,1.0],xtitle='<c>',ytitle='Alpha!I200!N',/xlog,/ylog
;plot,indgen(2),indgen(2),/nodata,xrange=[0.,0.12],yrange=[0.,0.2],ytitle='Alpha!I200!N',xtitle='<c>'
;for s=0,4 do begin 
;oplot,cmean(*),ratio(*,s),psym=sim(s+1),thick=2,col=color(s+1)
;endfor
;oplot,cmean(*),ratios(*),psym=sim(0),thick=2,col=color(0)
;for s=0,4 do begin 
;oplot,ratio(*,s),cmean(*),psym=4,thick=2,col=color(s+1)
;oplot,ratios(*),cmean(*),psym=4,thick=2,col=color(0)
;endfor
;legend,['R200/3','60','100','200','400','600'],psym=[sim(*)],col=[color(*)],/right_legend,/bottom_legend
;device,/close

;set_plot,'x'
;window,1
set_plot,'ps'
loadct,13
device,filename=fold2+'w-alpha_2.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[1e-3,0.2],yrange=[1e-2,1.0],xtitle='<w>',ytitle='Alpha!I200!N'
;plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[1e-3,0.2],yrange=[1e-2,1.0],xtitle='Alpha!I200!N',ytitle='<w>'
for s=0,4 do begin 
oplot,wmean(*),ratio(*,s),psym=sim(s+1),thick=2,col=color(s+1)
endfor
oplot,wmean(*),ratios(*),psym=sim(0),thick=2,col=color(0)
;for s=0,4 do begin 
;oplot,ratio(*,s),wmean(*),psym=4,thick=2,col=color(s+1)
;oplot,ratios(*),wmean(*),psym=4,thick=2,col=color(0)
;endfor
legend,['R200/3','60','100','200','400','600'],psym=[sim(*)],col=[color(*)],/right_legend,/bottom_legend
device,/close


end

;......STUDIO DELLA RELAZIONE ALPHA-SCALA DI FILTRAGGIO.....

pro scala_kolmo

fold2='/home/STUDENTI/matteo.angelinelli/Procedure/'
fold_out='/home/STUDENTI/matteo.angelinelli/Output/'
readcol,fold2+'kolmo.dat',sc,a0,err0,a2,err2,format='(f,f,f,f,f)'

nc=size(sc)

alfa=a0+a2
err=err0+err2
;err=stddev(alfa)

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'kolmogorov_2.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[50.,700.],yrange=[1e-3,1.],ytitle='Alpha!I200!N',xtitle='Smoothing scale [kpc]',/ystyle,/xstyle,/xlog,/ylog
oplot,sc,alfa,psym=4,thick=2,col=50
w=fltarr(nc(1))

w(*)=1.;/err(*)^2.

a=[8.e-4,0.8]

yfit=curvefit(sc,alfa,w,a,sigma,function_name='funct_kolmo',itmax=1000,chisq=chisq,yerror=yerr,iter=niter,status=status)

a1=[1.,-5./3.]

yfit2=curvefit(3.14*2./sc,alfa,w,a1,sigma1,function_name='funct_kolmo',itmax=1000,chisq=chisq1,yerror=yerr1,iter=niter1,status=status1)

print,a
print,sigma
print,chisq
print,status

chi=0.
for i=0,nc(1)-1 do begin 
chi+=(alfa(i)-yfit(i))^2./err(i)^2.
endfor

print,chi

print,a1
print,sigma1
print,chisq1
print,status1

chi1=0.
for i=0,nc(1)-1 do begin 
chi1+=(alfa(i)-yfit2(i))^2./err(i)^2.
endfor

print,chi1

;close,1
;openw,1,fold_out+'fit_parameters_smoothing.dat'
;printf,1,'To make fit we use the function: y=a(0)*x^a(1)' 
;print,status
;printf,1,'The fitting parameters are: ',a(0),a(1)
;printf,1,yfit
;printf,1,a0a2
;close,1

oplot,sc,yfit,linestyle=0,thick=2,col=250
;xyouts,100,0.13,'We use the function: alpha = a(0) * smoothing_scale!Ea(1)!N'
;xyouts,100,0.125,'We obtain for a(0)+-err(0)'
;xyouts,100,0.12,string(a(0),sigma(0))
;xyouts,100,0.115,'And for a(1)+-err(1):'
;xyouts,100,0.11,string(a(1),sigma(1))

legend,['Data','Fit'],linestyle=['',0],psym=[4,''],col=[50,250],/right_legend,/bottom_legend
device,/close

;window,1
set_plot,'ps'
loadct,13
device,filename=fold_out+'kolmogorov_3.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[5.e-3,0.15],yrange=[5e-4,1.5e-1],/xlog,/ylog,ytitle='Alpha(r!I200!N)',xtitle='k (=2*pi/l) [kpc!E-1!N]',/ystyle,/xstyle
oplot,2*3.14/sc,yfit2,linestyle=3,thick=2,col=250
oplot,2*3.14/sc,alfa,psym=4,thick=2,col=50
legend,['Data','Fit'],linestyle=['',0],psym=[4,''],col=[50,250],/right_legend,/bottom_legend

device,/close

end

pro funct_kolmo, x,a,f,pder      

    f=a[0]*x^a[1] 
    pder=[[x^a[1]],[a[0]*a[1]*x^(a[1]-1)]]
end

;.....CALCOLO DELLA MASSA TOTALE.......

function mass_tot,res,r200,overd
rhoc=10.^(-29.063736)
vol=3.*alog10(res)+3.*alog10(3.086*10.^21.)
m_tot=alog10((4./3.)*!pi)+3.*alog10(r200)+vol+alog10(rhoc*overd)

return,m_tot
end

;....CALCOLO DELL MASSA IDROSTATICA......

function mass_hydP,r,pma,pd,r200,res

r=r*r200

g=6.67*1e-8
kpc_cm=3.086*10.^(21.)
rad=r*kpc_cm;/float(res)   <---CHECK
r200=r200*kpc_cm;/float(res)  <---CHECK
nr=size(rad)
mass=fltarr(nr(1))
derp=deriv(pma)

for i=1,nr(1)-1 do begin
if rad(i) le r200 then begin 
   mass(i)=-1*(1e-20*(1e-20*rad(i))/(1e-20*kpc_cm))*rad(i)*derp(i)/float(G*pd(i)*res)
  ; print,"hydro test new",mass(i),rad(i),derp(i),G,pd(i)
   r200_index=i
endif
if rad(i) gt r200 then begin 
goto,finish
endif
endfor
finish:

mass200=20.+alog10(mass(r200_index))

return,mass200
end

pro alphamass

sc=['3','5','10','20','30']
fold_out='/home/STUDENTI/matteo.angelinelli/Procedure/'
fold2='/home/STUDENTI/matteo.angelinelli/Output/'

readcol,fold2+'catalog_cluster_OK.dat',id1,r200,mhyd,mtot,deltam,alfa,format='(a,f,f,f,f,f)',skipline=1

pos=where(deltam ge 0,nr)

ID=id1(pos)

n=size(ID)

ratio=fltarr(n(1),5)
ratios=fltarr(n(1))

for s=0,4 do begin 
for c=0,n(1)-1 do begin 
  readcol,fold2+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+ID(c)+'_pressure_30_'+sc(s)+'.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

  rv2=where(rr1 gt 0.99 and rr1 lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(rr1 gt 0.98 and rr1 lt 1.02,nr) 
  rv=rv2(0)                            
endif
 
  ratio(c,s)=rratio1(rv)

endfor
endfor

for c=0,n(1)-1 do begin 
  readcol,fold2+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+ID(c)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

  rv2=where(rr2 gt 0.99 and rr2 lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(rr2 gt 0.98 and rr2 lt 1.02,nr) 
  rv=rv2(0)                            
endif

  ratios(c)=rratio2(rv)

endfor

color=indgen(6)*55.
sim=[4,5,6,4,5,6]

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'dm-alpha_2.eps',/color
plot,deltam(pos),alfa(pos),psym=4,xtitle='(Mtot-Mhyd)/Mtot',ytitle='Alpha!I200!N',/xlog,/ylog
device,/close

;window,1
set_plot,'ps'
loadct,13
device,filename=fold_out+'alphaour-alphaeckert_2.eps',/color
plot,ratios,alfa(pos),psym=sim(0),col=color(0),xtitle='Alpha!I200,Our!N',ytitle='Alpha!I200,Eckert!N',/xlog,/ylog,xrange=[1e-3,1]
for s=0,4 do begin 
oplot,ratio(*,s),alfa(pos),psym=sim(s+1),col=color(s+1)
endfor
legend,['R200/3','60','100','200','400','600'],psym=[sim(*)],col=[color(*)],/right_legend,/bottom_legend
device,/close

;window,4
set_plot,'ps'
loadct,13
device,filename=fold_out+'dm-alphaour_2.eps',/color
plot,deltam(pos),ratios,psym=sim(0),col=color(0),ytitle='Alpha!I200,Our!N',xtitle='(Mtot-Mhyd)/Mtot',/xlog,/ylog,yrange=[1e-3,1]
for s=0,4 do begin 
oplot,deltam(pos),ratio(*,s),psym=sim(s+1),col=color(s+1)
endfor
legend,['R200/3','60','100','200','400','600'],psym=[sim(*)],col=[color(*)],/right_legend,/bottom_legend
device,/close

end

pro plot_eckert

fold_out='/home/STUDENTI/matteo.angelinelli/Output/'
fold_file='/home/STUDENTI/matteo.angelinelli/Procedure/'

cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]
sc=['3','5','10','20','30']
nc1=size(cl_list)
nbr=30

r_plot1=fltarr(nbr,nc1(1))
ratio_plot1=fltarr(nbr,nc1(1))
r_plot=fltarr(nbr)
ratio_plot=fltarr(nbr)
r_plt1=fltarr(nbr,nc1(1))
ratio_plt1=fltarr(nbr,nc1(1))
r_plt=fltarr(nbr)
ratio_plt=fltarr(nbr)
color=indgen(6)*55.
line=indgen(6)
line(*)=2

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'Plot_median_eckert.eps',/color
plot,indgen(2),indgen(2),/nodata,xrange=[-0.1,2.0],yrange=[1e-4,1.],xtitle='r/r!d200!n',ytitle='Alpha',/xstyle,/ystyle,/ylog

for c=0,nc1(1)-1 do begin 
  readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(c)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

r_plt1(*,c)=rr2
ratio_plt1(*,c)=rratio2
endfor

for i=0,nbr-1 do begin 
r_plt(i)=median(r_plt1(i,*))
ratio_plt(i)=median(ratio_plt1(i,*))
endfor

oplot,r_plt,ratio_plt,linestyle=2,thick=2,col=color(0)


for s=0,4 do begin 
for c=0,nc1(1)-1 do begin 
  ;readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(c)+'_pressure_scalamobile.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

readcol,fold_out+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(c)+'_pressure_30_'+sc(s)+'.dat',rr,rratio,format='(f,f)',/silent ;...file with single cluster

r_plot1(*,c)=rr
ratio_plot1(*,c)=rratio
endfor

for i=0,nbr-1 do begin 
r_plot(i)=median(r_plot1(i,*))
ratio_plot(i)=median(ratio_plot1(i,*))
endfor

oplot,r_plot,ratio_plot,linestyle=2,thick=2,col=color(s+1)

endfor

readcol,fold_file+'eckert_red_error_K0.dat',id,m200,alfa200,alfa200p,alfa200m,alfa500,alfa500p,alfa500m,k0,format='(A,F,F,F,F,F,F,F,F)',skipline=2,/silent

alfa200=alfa200/100.
alfa200p=alfa200p/100.
alfa200m=alfa200m/100.
alfa500=alfa500/100.
alfa500p=alfa500p/100.
alfa500m=alfa500m/100.

nec=size(id)
r200plot=fltarr(nec(1))
r500plot=fltarr(nec(1))
r200plot(*)=0.99+4e-2*indgen(nec(1))/float(nec(1))
r500plot(*)=0.69+4e-2*indgen(nec(1))/float(nec(1))

oploterr,r200plot,alfa200,alfa200p,5
oploterr,r200plot,alfa200,alfa200m,5
oploterr,r500plot,alfa500,alfa500p,4
oploterr,r500plot,alfa500,alfa500m,4

legend,['Eckert et al. +18 -- R500','Eckert et al. +18 -- R200'],psym=[4,5],/left,/bottom
legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],linestyle=line,col=color,/right_legend,/bottom_legend

device,/close

end


pro profile

fold_out='/home/STUDENTI/matteo.angelinelli/Output/'

cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]

nc=size(cl_list)


colorcluster=indgen(nc(1))*20.
line=indgen(nc(1))
line=0.

set_plot,'ps'
loadct,13
device,filename=fold_out+'Temperature_profile.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[1.e-3,12.],yrange=[1.e4,1.e8],xtitle='r/r!d200!N',ytitle='Temperature [K]' 

for c=0,nc(1)-1 do begin
readcol,fold_out+'profile_d_dm_t_'+cl_list(c)+'_OK.dat',r,d,dm,t,format='(f,f,f,f)',skipline=1

oplot,r,t,col=colorcluster(c)

endfor

legend,[cl_list],linestyle=line,col=colorcluster,/left_legend,/bottom_legend
device,/close

end


pro mappe
  foldc='/home/data/DATA/ISC/'
fold_out='/home/STUDENTI/matteo.angelinelli/Output/'
cluster='IT90_3'
snapc='193'

 CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine
 
  fold=foldc+cluster+'/'
  snap='0'+string(snapc,'(i3)')
  file1=fold+'/declust_z'+snap   ;...the field names are Density, Temperature, Dark_Matter_Densiy
  file2=fold+'declust_v_z'+snap ;..." x-velocity, y-velocity, z-velocity

  if cluster eq 'IT3' or cluster eq 'IT7' or cluster eq 'IT6' or cluster eq 'IT15B' then file2=file1  ;...these files where written in a slightly different way

  fileconv='deDD'+snap+'.conv2' ;...conversion factors from enzo's internal units to cgs

  openr,3,fold+fileconv
  readf,3,snapshot_id
  readf,3,time
  readf,3,redshift
  readf,3,convd  ;...conversion factor from enzo's to g/cm^3
  readf,3,convv   ;...      "      "        cm/s
  close,3

  ;....the conversion factors change with redshift and simulation, must be found into enzo output files or are written elsewhere during the data reconstruction (ask Franco)
 ;....HDF5 (see if it's easy to download the hdf5 libraries for shell commands
  ;...https://www.hdfgroup.org/downloads/


  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Density')
  datasp=h5d_get_space(itemp)
  d = H5D_READ(itemp)
  h5d_close,itemp

  
  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Dark_Matter_Density')
  datasp=h5d_get_space(itemp)
    dm = H5D_READ(itemp)
  h5d_close,itemp

  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Temperature')
  datasp=h5d_get_space(itemp)
  temp = H5D_READ(itemp)
  h5d_close,itemp

  ido=h5f_open(file2)
  itemp=h5d_open(ido,'x-velocity')
  datasp=h5d_get_space(itemp)
  vx = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'y-velocity')
  datasp=h5d_get_space(itemp)
  vy = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'z-velocity')
  datasp=h5d_get_space(itemp)
  vz = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido

d=d*convd
dm=dm*convd

;....some manipulation of grid size needed only in the rare case the
;HDF5 dataset is not a cube

  n3=size(d)                    ;...this gives the size of the density field just read
  print,n3

  n0=n3(1)
  n1=n3(2)
  n2=n3(3)
   n=n0
  if n0 ne n1 or n0 ne n2 or n1 ne n2 then begin
     print,"making the dataset cubic"
                                ;....probably there is a much more
                                ;elegant way of doing the following, anyways....
     n=min(n3(1:3))
     dnew=d(0:n-1,0:n-1,0:n-1)
     tnew=temp(0:n-1,0:n-1,0:n-1)
     dmnew=dm(0:n-1,0:n-1,0:n-1)
     vxnew=vx(0:n-1,0:n-1,0:n-1)
     vynew=vy(0:n-1,0:n-1,0:n-1)
     vznew=vz(0:n-1,0:n-1,0:n-1)
     d=dnew
     dm=dmnew
     temp=tnew
     vx=vxnew
     vy=vynew
     vz=vznew
     vxnew=0
     vynew=0
     vznew=0
  end   

sc=10.

  vx_smooth=smooth(vx,sc)
  vy_smooth=smooth(vy,sc)
  vz_smooth=smooth(vz,sc)

  vturbx=vx-vx_smooth
  vturby=vy-vy_smooth
  vturbz=vz-vz_smooth

  vturb=convv*sqrt(vturbx^2.+vturby^2.+vturbz^2.)  ;...velocity module in cm/s
  vtot=convv*sqrt(vx^2.+vy^2.+vz^2.)
  vsmooth=convv*sqrt(vx_smooth^2.+vy_smooth^2.+vz_smooth^2.)  


vx_smooth=0
vy_smooth=0
vz_smooth=0
vturbx=0
vturby=0
vturbz=0

vol=3.*alog10(20.)+3.*alog10(3.086*10.^21.)
;...the first 3 are the input un-filtered velocities, assumed to be in cgs
;...the input density must be in cgs
;...vol=log10(cell volume in cgs)
;...resdhift
;....flux is the (output) kinetic energy flux through shocked cells, in log10(erg/s)
;....macx is the (output) Mach number of shocked cells
;....div is the (output) 3D divergence of the velocity field, in cm/s/cell

a=shock_finder(convv*vx,convv*vy,convv*vz,d,temp,vol-3.*alog10(1.+redshift),flux,macx,div) 



n=320.
vect=fltarr(n,n,2)
v_vel=fltarr(n,n,3)

  for k=(n/2.-1)-20.,(n/2.-1)+20. do begin ;...to make the vector which contains the projections along the z-axis of the various velocities 
  vect(*,*,0)+=vturb(*,*,k)
  endfor

im=where(macx gt 1.3,nr)

vturb(im)=0.

  for k=(n/2.-1)-20.,(n/2.-1)+20. do begin ;...to make the vector which contains the projections along the z-axis of the various velocities 
  vect(*,*,1)+=vturb(*,*,k)
  endfor

vect(*,*,*)=vect(*,*,*)/40.

writefits,fold_out+'maps_turbmach_central.fits',alog10(vect(*,*,*)) 

 for k=0,n-1 do begin 
 v_vel(*,*,0)+=vtot(*,*,k)
 v_vel(*,*,1)+=vsmooth(*,*,k)
 v_vel(*,*,2)+=vturb(*,*,k)
endfor

v_vel(*,*,*)/=n

writefits,fold_out+'maps_vel_total.fits',alog10(v_vel(*,*,*)) 

v_mach=fltarr(n,n,2)

for k=(n/2.-1)-10.,(n/2.-1)+10. do begin 
v_mach(*,*,0)+=macx(*,*,k)
v_mach(*,*,1)+=10.^(flux(*,*,k)-10.)
endfor

v_mach(*,*,*)/=20.
v_mach(*,*,1)=alog10(v_mach(*,*,1))+10.

writefits,fold_out+'maps_mach_slice_2.fits',v_mach(*,*,*)

end

pro profili_vturb

  foldc='/home/data/DATA/ISC/'
fold_out='/home/STUDENTI/matteo.angelinelli/Output/'

 CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine

cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]
sn_list=[199,196,195,193,199,243,228,242,102,100,115,167,106,109] 

nc=size(cl_list)

colorcluster=indgen(nc(1))*20.
line=indgen(nc(1))
line=0.

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'Profile_velocity_turb.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[1.e-3,1.5e1],yrange=[1.e6,1.e8],xtitle='r/r!d200!N',ytitle='Velocity [cm s!E-1!N]' 

pveltot=fltarr(280,nc(1))
pvelftot=fltarr(280,nc(1))
pvelnftot=fltarr(280,nc(1))

for c=0,nc(1)-1 do begin
;c=0
cluster=cl_list(c)
snapc=sn_list(c)

if c eq 11 then begin 
res=10.
endif

if c ne 11 then begin 
res=20.
endif

  fold=foldc+cluster+'/'
  snap='0'+string(snapc,'(i3)')
  file1=fold+'/declust_z'+snap   ;...the field names are Density, Temperature, Dark_Matter_Densiy
  file2=fold+'declust_v_z'+snap ;..." x-velocity, y-velocity, z-velocity

  if cluster eq 'IT3' or cluster eq 'IT7' or cluster eq 'IT6' or cluster eq 'IT15B' then file2=file1  ;...these files where written in a slightly different way

  fileconv='deDD'+snap+'.conv2' ;...conversion factors from enzo's internal units to cgs

  openr,3,fold+fileconv
  readf,3,snapshot_id
  readf,3,time
  readf,3,redshift
  readf,3,convd  ;...conversion factor from enzo's to g/cm^3
  readf,3,convv   ;...      "      "        cm/s
  close,3

  ;....the conversion factors change with redshift and simulation, must be found into enzo output files or are written elsewhere during the data reconstruction (ask Franco)
 ;....HDF5 (see if it's easy to download the hdf5 libraries for shell commands
  ;...https://www.hdfgroup.org/downloads/


  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Density')
  datasp=h5d_get_space(itemp)
  d = H5D_READ(itemp)
  h5d_close,itemp

  
  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Dark_Matter_Density')
  datasp=h5d_get_space(itemp)
    dm = H5D_READ(itemp)
  h5d_close,itemp

  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Temperature')
  datasp=h5d_get_space(itemp)
  temp = H5D_READ(itemp)
  h5d_close,itemp

  ido=h5f_open(file2)
  itemp=h5d_open(ido,'x-velocity')
  datasp=h5d_get_space(itemp)
  vx = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'y-velocity')
  datasp=h5d_get_space(itemp)
  vy = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'z-velocity')
  datasp=h5d_get_space(itemp)
  vz = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido

;....some manipulation of grid size needed only in the rare case the
;HDF5 dataset is not a cube

  n3=size(d)                    ;...this gives the size of the density field just read
  print,n3

  n0=n3(1)
  n1=n3(2)
  n2=n3(3)
   n=n0
  if n0 ne n1 or n0 ne n2 or n1 ne n2 then begin
     print,"making the dataset cubic"
                                ;....probably there is a much more
                                ;elegant way of doing the following, anyways....
     n=min(n3(1:3))
     dnew=d(0:n-1,0:n-1,0:n-1)
     tnew=temp(0:n-1,0:n-1,0:n-1)
     dmnew=dm(0:n-1,0:n-1,0:n-1)
     vxnew=vx(0:n-1,0:n-1,0:n-1)
     vynew=vy(0:n-1,0:n-1,0:n-1)
     vznew=vz(0:n-1,0:n-1,0:n-1)
     d=dnew
     dm=dmnew
     temp=tnew
     vx=vxnew
     vy=vynew
     vz=vznew
     vxnew=0
     vynew=0
     vznew=0
  end  

d=d*convd
dm=dm*convd 

sc=10.

  vx_smooth=smooth(vx,sc)
  vy_smooth=smooth(vy,sc)
  vz_smooth=smooth(vz,sc)

  vturbx=vx-vx_smooth
  vturby=vy-vy_smooth
  vturbz=vz-vz_smooth

  vturb=convv*sqrt(vturbx^2.+vturby^2.+vturbz^2.)  ;...velocity module in cm/s
  vtot=convv*sqrt(vx^2.+vy^2.+vz^2.)
  vsmooth=convv*sqrt(vx_smooth^2.+vy_smooth^2.+vz_smooth^2.)  


vx_smooth=0
vy_smooth=0
vz_smooth=0
vturbx=0
vturby=0
vturbz=0

vol=3.*alog10(res)+3.*alog10(3.086*10.^21.)
;...the first 3 are the input un-filtered velocities, assumed to be in cgs
;...the input density must be in cgs
;...vol=log10(cell volume in cgs)
;...resdhift
;....flux is the (output) kinetic energy flux through shocked cells, in log10(erg/s)
;....macx is the (output) Mach number of shocked cells
;....div is the (output) 3D divergence of the velocity field, in cm/s/cell

a=shock_finder(convv*vx,convv*vy,convv*vz,d,temp,vol-3.*alog10(1.+redshift),flux,macx,div) 

 xm=where(d+dm eq max(d+dm))
  xc=array_indices(d,xm)

  print,"the center of this cluster is",xc

  rx=indgen(n)-xc(0)
  ry=indgen(n)-xc(1)
  rz=indgen(n)-xc(2)
  rx=reform(rx,n,1,1)
  ry=reform(ry,1,n,1)
  rz=reform(rz,1,1,n)
  ra=rebin(rx,n,n,n)^2
  ra=temporary(ra)+rebin(ry,n,n,n)^2                 ; conserve memory
  ra=temporary(ra)+rebin(rz,n,n,n)^2
  ra=fix(temporary(sqrt(ra)))  ;...array of 3D radial distances

pvel=fltarr(max(ra))
pvelf=fltarr(max(ra))
pvelnf=fltarr(max(ra))
pd=fltarr(max(ra))
pdm=fltarr(max(ra))

for r=0,max(ra)-1 do begin
  
   rw=where(ra ge r and ra lt r+1,nw)   ;…we can also do r
   rk=where(ra ge r and ra lt r+1 and macx lt 1.3,nk) 
   rz=where(ra ge r and ra lt r+1 and macx gt 1.3,nz)

   pvel(r)=mean(vturb(rw))
   pvelf(r)=mean(vturb(rk))
   pvelnf(r)=mean(vturb(rz))

    pd(r)=mean(d(rw))
    pdm(r)=mean(dm(rw))
    
endfor

x=res*(0.5+indgen(max(ra)))

overd=200 ;....overdensity 
r200=rvir(x,pd,pdm,overd)   ;...procedure to derive the radius enclosing a >overd overdensity
x=x/float(r200*res)

oplot,x,pvel,linestyle=0,col=colorcluster(c)
oplot,x,pvelf,linestyle=2,col=colorcluster(c)
oplot,x,pvelnf,linestyle=3,col=colorcluster(c)

pveltot(*,c)=pvel
pvelftot(*,c)=pvelf
pvelnftot(*,c)=pvelnf
endfor

legend,[cl_list],linestyle=line,col=colorcluster,/left_legend,/bottom_legend
legend,['All Cells', 'No Shocked Cells','Shocked Cells'],linestyle=[0,2,3],/center,/bottom_legend
device,/close 

pvelm=fltarr(280)
pvelfm=fltarr(280)
pvelnfm=fltarr(280)

for i=0,279 do begin 
pvelm(i)=mean(pveltot(i,*))
pvelfm(i)=mean(pvelftot(i,*))
pvelnfm(i)=mean(pvelnftot(i,*))
endfor

set_plot,'ps'
loadct,13
device,filename=fold_out+'Profile_velocity_turb_median.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[1.e-3,1.5e1],yrange=[1.e6,1.e8],xtitle='r/r!d200!N',ytitle='Velocity [cm s!E-1!N]'
oplot,x,pvelm,linestyle=0,thick=2
oplot,x,pvelfm,linestyle=2,thick=2
oplot,x,pvelnfm,linestyle=3,thick=2

legend,['Median All Cells','Median No Shocked Cells','Median Shocked Cell'],linestyle=[0,2,3],/right,/bottom
device,/close
end

pro profili_cumulativi

t0=systime(1)

  foldc='/home/data/DATA/ISC/'
  fold_out='/home/STUDENTI/matteo.angelinelli/Output/'    ;...output folder

cluster='IT90_3'
snapc='193'

  CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine
 
  fold=foldc+cluster+'/'
  snap='0'+string(snapc,'(i3)')
  file1=fold+'/declust_z'+snap   ;...the field names are Density, Temperature, Dark_Matter_Densiy
  file2=fold+'declust_v_z'+snap ;..." x-velocity, y-velocity, z-velocity

  if cluster eq 'IT3' or cluster eq 'IT7' or cluster eq 'IT6' or cluster eq 'IT15B' then file2=file1  ;...these files where written in a slightly different way

  fileconv='deDD'+snap+'.conv2' ;...conversion factors from enzo's internal units to cgs

  openr,3,fold+fileconv
  readf,3,snapshot_id
  readf,3,time
  readf,3,redshift
  readf,3,convd  ;...conversion factor from enzo's to g/cm^3
  readf,3,convv   ;...      "      "        cm/s
  close,3

  print,cluster,redshift,convd,convv

  ;....the conversion factors change with redshift and simulation, must be found into enzo output files or are written elsewhere during the data reconstruction (ask Franco)
 ;....HDF5 (see if it's easy to download the hdf5 libraries for shell commands
  ;...https://www.hdfgroup.org/downloads/


  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Density')
  datasp=h5d_get_space(itemp)
  dens= H5D_READ(itemp)
  h5d_close,itemp

  
  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Dark_Matter_Density')
  datasp=h5d_get_space(itemp)
  dens_dm = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'Temperature')
  datasp=h5d_get_space(itemp)
  temp = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido


  ido=h5f_open(file2)
  itemp=h5d_open(ido,'x-velocity')
  datasp=h5d_get_space(itemp)
  vx = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'y-velocity')
  datasp=h5d_get_space(itemp)
  vy = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'z-velocity')
  datasp=h5d_get_space(itemp)
  vz = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido

dens=convd*dens       ;...density in g/cm^3 
dens_dm=convd*dens_dm

;....some manipulation of grid size needed only in the rare case the
;HDF5 dataset is not a cube

  n3=size(dens)                    ;...this gives the size of the density field just read
  print,n3

  n0=n3(1)
  n1=n3(2)
  n2=n3(3)
   n=n0
  if n0 ne n1 or n0 ne n2 or n1 ne n2 then begin
     print,"making the dataset cubic"
                                ;....probably there is a much more
                                ;elegant way of doing the following, anyways....
     n=min(n3(1:3))
     dnew=dens(0:n-1,0:n-1,0:n-1)
     tnew=temp(0:n-1,0:n-1,0:n-1)
     dmnew=dens_dm(0:n-1,0:n-1,0:n-1)
     vxnew=vx(0:n-1,0:n-1,0:n-1)
     vynew=vy(0:n-1,0:n-1,0:n-1)
     vznew=vz(0:n-1,0:n-1,0:n-1)
     dens=dnew
     dens_dm=dmnew
     temp=tnew
     vx=vxnew
     vy=vynew
     vz=vznew
     vxnew=0
     vynew=0
     vznew=0
  end   

t1=systime(1)
print, '...read complete in',(systime(1)-t0),' seconds...'

vol=3.*alog10(20.)+3.*alog10(3.086*10.^21.)
;...the first 3 are the input un-filtered velocities, assumed to be in cgs
;...the input density must be in cgs
;...vol=log10(cell volume in cgs)
;...resdhift
;....flux is the (output) kinetic energy flux through shocked cells, in log10(erg/s)
;....macx is the (output) Mach number of shocked cells
;....div is the (output) 3D divergence of the velocity field, in cm/s/cell

;a=shock_finder(convv*vx,convv*vy,convv*vz,dens,temp,vol-3.*alog10(1.+redshift),flux,macx,div) 

t2=systime(1)
print, '...shocks finder complete in',(systime(1)-t1),' seconds...'
  
  xm=where(dens+dens_dm eq max(dens+dens_dm))
  xc=array_indices(dens,xm)

  print,"the center of this cluster is",xc

  rx=indgen(n)-xc(0)
  ry=indgen(n)-xc(1)
  rz=indgen(n)-xc(2)
  rx=reform(rx,n,1,1)
  ry=reform(ry,1,n,1)
  rz=reform(rz,1,1,n)
  ra=rebin(rx,n,n,n)^2
  ra=temporary(ra)+rebin(ry,n,n,n)^2                 ; conserve memory
  ra=temporary(ra)+rebin(rz,n,n,n)^2
  ra=fix(temporary(sqrt(ra)))  ;...array of 3D radial distances


print, '...grid complete in ',systime(1)-t2,' seconds...'
t3=systime(1)

n=320.


sc=30
print, sc

vtot=fltarr(n,n,n)
vturb=fltarr(n,n,n)

  vx_smooth=smooth(vx,sc)
  vy_smooth=smooth(vy,sc)
  vz_smooth=smooth(vz,sc)

  vturbx=vx-vx_smooth
  vturby=vy-vy_smooth
  vturbz=vz-vz_smooth

  vmod_turb=convv*sqrt(vturbx^2.+vturby^2.+vturbz^2.)  ;...velocity module in cm/s
  vmod_tot=convv*sqrt(vx^2.+vy^2.+vz^2.)
vx_smooth=0
vy_smooth=0
vz_smooth=0
vturbx=0
vturby=0
vturbz=0

t4=systime(1)
print, '...velocity fields complete in',systime(1)-t3,' seconds...'

  ekt=fltarr(max(ra))
  eth=fltarr(max(ra))
  ektot=fltarr(max(ra))
  kb=1.38*10.^(-16.)
  mh=1.67*10.^(-24.)
  mu=0.59

  prof=fltarr(max(ra),3)
  pd=fltarr(max(ra))
  pdm=fltarr(max(ra))
   
  for r=0,max(ra)-1 do begin
  
    rk=where(ra ge r and ra lt r+1,nm)   ;…we can also do rk=where(ra ge r and ra lt r+2,nr) for example
    
    ekt(rk)=0.5*dens(rk)*vmod_turb(rk)^(2.) 
    eth(rk)=(3./2.)*kb*temp(rk)*dens(rk)/(mh*mu) 
    ektot(rk)=0.5*dens(rk)*vmod_tot(rk)^(2.)

    prof(r,0)=total(ekt(rk))        ;....cumulative energy in shell
    prof(r,1)=total(eth(rk))
    prof(r,2)=total(ektot(rk))
    pd(r)=mean(dens(rk))
    pdm(r)=mean(dens_dm(rk))
    
 endfor
 
   ekt=0
   eth=0  
   ektot=0

   prof(*,*)=alog10(prof(*,*))+vol

 pcumul=fltarr(max(ra),3)   
   
  for j=0,2 do begin  
   for i=1,max(ra)-1 do begin 
  pcumul(i,j)=total(10.^(prof(0:i-1,j)-30))    ;...cumulative energy from 0 to i-1
   endfor
  endfor 

prof=0

x=20.*(0.5+indgen(max(ra)))

overd=200 ;....overdensity 
res=20 ;...kpc
r200=rvir(x,pd,pdm,overd)   ;...procedure to derive the radius enclosing a >overd overdensity
x=x/float(r200*res)

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'Cumulative_profile_30.eps',/color
plot,indgen(2),indgen(2),/nodata,yrange=[54.,62.],xrange=[5.e-2,1.e1],/xlog,col=0,xtitle='r/r!d200!n',ytitle='Energy [log!d10!n(erg)] '
oplot,x,alog10(pcumul(*,0))+30.,col=0
oplot,x,alog10(pcumul(*,1))+30.,col=100
oplot,x,alog10(pcumul(*,2))+30.,col=150
legend,['Thermal Energy','Kinetic Energy','Turbulent Energy'],linestyle=[0,0,0],col=[100,150,0],/right_legend,/bottom_legend
legend,['600 kpc'],col=[0],/left_legend,/bottom_legend
device,/close

end


pro isto

 foldc='/home/data/DATA/ISC/'
  fold_out='/home/STUDENTI/matteo.angelinelli/Output/'    ;...output folder

cluster='IT90_3'
snapc='193'

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold_out+'Mach_isto.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xtitle='Mach Number',ytitle='Numbers of Cells/Binsize',xrange=[1,1e3],yrange=[1.e4,1.e10]


  CPU,TPOOL_NTHREADS=32 ;.....to accelerate some calculations on this machine
 
  fold=foldc+cluster+'/'
  snap='0'+string(snapc,'(i3)')
  file1=fold+'/declust_z'+snap   ;...the field names are Density, Temperature, Dark_Matter_Densiy
  file2=fold+'declust_v_z'+snap ;..." x-velocity, y-velocity, z-velocity

  if cluster eq 'IT3' or cluster eq 'IT7' or cluster eq 'IT6' or cluster eq 'IT15B' then file2=file1  ;...these files where written in a slightly different way

  fileconv='deDD'+snap+'.conv2' ;...conversion factors from enzo's internal units to cgs

  openr,3,fold+fileconv
  readf,3,snapshot_id
  readf,3,time
  readf,3,redshift
  readf,3,convd  ;...conversion factor from enzo's to g/cm^3
  readf,3,convv   ;...      "      "        cm/s
  close,3

  print,cluster,redshift,convd,convv

  ;....the conversion factors change with redshift and simulation, must be found into enzo output files or are written elsewhere during the data reconstruction (ask Franco)
 ;....HDF5 (see if it's easy to download the hdf5 libraries for shell commands
  ;...https://www.hdfgroup.org/downloads/


  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Density')
  datasp=h5d_get_space(itemp)
  dens= H5D_READ(itemp)
  h5d_close,itemp

  
  ido=h5f_open(file1)
  itemp=h5d_open(ido,'Dark_Matter_Density')
  datasp=h5d_get_space(itemp)
  dens_dm = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'Temperature')
  datasp=h5d_get_space(itemp)
  temp = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido


  ido=h5f_open(file2)
  itemp=h5d_open(ido,'x-velocity')
  datasp=h5d_get_space(itemp)
  vx = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'y-velocity')
  datasp=h5d_get_space(itemp)
  vy = H5D_READ(itemp)
  h5d_close,itemp

  itemp=h5d_open(ido,'z-velocity')
  datasp=h5d_get_space(itemp)
  vz = H5D_READ(itemp)
  h5d_close,itemp
  h5f_close,ido

dens=convd*dens       ;...density in g/cm^3 
dens_dm=convd*dens_dm

;....some manipulation of grid size needed only in the rare case the
;HDF5 dataset is not a cube

  n3=size(dens)                    ;...this gives the size of the density field just read
  print,n3

  n0=n3(1)
  n1=n3(2)
  n2=n3(3)
   n=n0
  if n0 ne n1 or n0 ne n2 or n1 ne n2 then begin
     print,"making the dataset cubic"
                                ;....probably there is a much more
                                ;elegant way of doing the following, anyways....
     n=min(n3(1:3))
     dnew=dens(0:n-1,0:n-1,0:n-1)
     tnew=temp(0:n-1,0:n-1,0:n-1)
     dmnew=dens_dm(0:n-1,0:n-1,0:n-1)
     vxnew=vx(0:n-1,0:n-1,0:n-1)
     vynew=vy(0:n-1,0:n-1,0:n-1)
     vznew=vz(0:n-1,0:n-1,0:n-1)
     dens=dnew
     dens_dm=dmnew
     temp=tnew
     vx=vxnew
     vy=vynew
     vz=vznew
     vxnew=0
     vynew=0
     vznew=0
  end   


vol=3.*alog10(20.)+3.*alog10(3.086*10.^21.)
;...the first 3 are the input un-filtered velocities, assumed to be in cgs
;...the input density must be in cgs
;...vol=log10(cell volume in cgs)
;...resdhift
;....flux is the (output) kinetic energy flux through shocked cells, in log10(erg/s)
;....macx is the (output) Mach number of shocked cells
;....div is the (output) 3D divergence of the velocity field, in cm/s/cell

a=shock_finder(convv*vx,convv*vy,convv*vz,dens,temp,vol-3.*alog10(1.+redshift),flux,macx,div) 

 min_mach=1.3
 mim=alog10(min_mach)   
 mam=alog10(1100.)      
 ntm=100.  ;..number of mach bins
 bt=(mam-mim)/float(ntm)
 xm=10.^(mim+bt*indgen(ntm))

 imac=where(macx gt min_mach,n1)
 
hm=histogram(alog10(macx(imac)),min=mim,max=mam,binsize=bt)  
hm=hm/bt   


oplot,xm,hm/bt
device,/close

end

pro alfacw

fold_pro='/home/STUDENTI/matteo.angelinelli/Procedure/' 
fold2='/home/STUDENTI/matteo.angelinelli/Output/'

readcol,fold_pro+'parametriX.dat',id1,proi,c,w,format='(a,f,f,f)'

cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]
sc=['3','5','10','20','30']
nc=size(cl_list)

cmean=fltarr(nc(1))
wmean=fltarr(nc(1))


for k=0,nc(1)-1 do begin 

posid=where(id1 eq cl_list(k)) 

cmean(k)=mean(c(posid))
wmean(k)=mean(w(posid))

endfor

ratio=fltarr(nc(1),5)
ratios=fltarr(nc(1))
ratio0=fltarr(nc(1),5)
ratios0=fltarr(nc(1))


for s=0,4 do begin 
for k=0,nc(1)-1 do begin 
  readcol,fold2+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(k)+'_pressure_30_'+sc(s)+'.dat',rr1,rratio1,format='(f,f)',/silent ;...file with single cluster

  rv2=where(rr1 gt 0.99 and rr1 lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(rr1 gt 0.98 and rr1 lt 1.02,nr) 
  rv=rv2(0)                            
endif

  ratio(k,s)=rratio1(rv)
  ratio0(k,s)=rratio1(0)

endfor
endfor

for k=0,nc(1)-1 do begin 
  readcol,fold2+'Prof/Prof_inter/'+'Profile_interpolate_cluster_'+cl_list(k)+'_pressure_scalamobile.dat',rr2,rratio2,format='(f,f)',/silent ;...file with single cluster

  rv2=where(rr2 gt 0.99 and rr2 lt 1.01,nr) 
  rv=rv2(0)                            
if rv eq -1 then begin   
  rv2=where(rr2 gt 0.98 and rr2 lt 1.02,nr) 
  rv=rv2(0)                            
endif

  ratios(k)=rratio2(rv)
  ratios0(k)=rratio2(0)

endfor

sim=[4,5,6,4,5,6]
color=indgen(6)*55.

set_plot,'x'
window,0
plot,indgen(2),indgen(2),/nodata,xtitle='<c>',ytitle='alpha!I200!N',xrange=[1e-4,1],yrange=[1e-3,1],/xlog,/ylog
for s=0,4 do begin 
oplot,ratio0(*,s),cmean(*),psym=sim(s+1)
endfor
oplot,ratios0,cmean(*),psym=sim(0)
;set_plot,'x'
;window,0
;set_plot,'ps'
;loadct,13
;device,filename=fold2+'alpha-c-200.eps',/color
;plot,indgen(2),indgen(2),/nodata,xtitle='<c>',ytitle='alpha!I200!N',xrange=[5.e-2,0.2],yrange=[0.001,0.2],/xlog,/ylog
;for s=0,4 do begin 
;oplot,cmean(*),ratio(*,s),psym=sim(s+1),col=color(s+1)
;endfor
;oplot,cmean(*),ratios(*),psym=sim(0),col=color(0)
;legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],psym=sim,col=color,/right_legend,/bottom_legend
;device,/close

;window,1
;set_plot,'ps'
;loadct,13
;device,filename=fold2+'alpha-c-0.eps',/color
;plot,indgen(2),indgen(2),/nodata,xtitle='<c>',ytitle='alpha!I0!N',xrange=[5.e-2,0.2],yrange=[5.e-4,1.],/xlog,/ylog
;for s=0,4 do begin 
;oplot,cmean(*),ratio0(*,s),psym=sim(s+1),col=color(s+1)
;endfor
;oplot,cmean(*),ratios0(*),psym=sim(0),col=color(0)
;legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],psym=sim,col=color,/right_legend,/bottom_legend
;device,/close

;window,2
;set_plot,'ps'
;loadct,13
;device,filename=fold2+'alpha-w-200.eps',/color
;plot,indgen(2),indgen(2),/nodata,xtitle='<w>',ytitle='alpha!I200!N',xrange=[0.001,0.02],yrange=[0.001,0.2],/xlog,/ylog
;for s=0,4 do begin 
;oplot,wmean(*),ratio(*,s),psym=sim(s+1),col=color(s+1)
;endfor
;oplot,wmean(*),ratios(*),psym=sim(0),col=color(0)
;legend,['R200/3','60kpc','100kpc','200kpc','400kpc','600kpc'],psym=sim,col=color,/right_legend,/bottom_legend
;device,/close

end

pro massaparax

fold_pro='/home/STUDENTI/matteo.angelinelli/Procedure/' 
fold2='/home/STUDENTI/matteo.angelinelli/Output/'

readcol,fold_pro+'parametriX.dat',id1,proi,c,w,format='(a,f,f,f)'

cl_list=["IT90_0","IT90_1","IT90_2","IT90_3","IT90_4","IT92_0","IT92_1","IT92_2","IT1","IT3","IT7","IT10","IT62","IT6"]
sc=['3','5','10','20','30']
nc=size(cl_list)

cmean=fltarr(nc(1))
wmean=fltarr(nc(1))


for k=0,nc(1)-1 do begin 

posid=where(id1 eq cl_list(k)) 

cmean(k)=mean(c(posid))
wmean(k)=mean(w(posid))

endfor

readcol,fold2+'catalog_cluster_OK.dat',id,r200,mhyd,mtot,deltam,alpha200,format='(a,f,f,f,f,f)'

colorcluster=indgen(nc(1))*20.

;set_plot,'x'
;window,0
set_plot,'ps'
loadct,13
device,filename=fold2+'massac.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[5.e13,1.e15],yrange=[5.e-2,1.e-1],xtitle='Total Mass [M!Io!N]',ytitle='<c>'
plots,[mtot],[cmean],psym=4,col=[colorcluster]
legend,[cl_list],psym=4,col=[colorcluster],/right,/top
device,/close
;window,1
set_plot,'ps'
loadct,13
device,filename=fold2+'massaw.eps',/color
plot,indgen(2),indgen(2),/nodata,/xlog,/ylog,xrange=[5.e13,1.e15],yrange=[1.e-3,5.e-2],xtitle='Total Mass [M!Io!N]',ytitle='<w>'
plots,[mtot],[wmean],psym=4,col=[colorcluster]
legend,[cl_list],psym=4,col=[colorcluster],/right,/top
device,/close
end

;+

FUNCTION POLY2D_PROJECT,IMAGE,MASK,ORDER,BASIS=BASIS,GIJ=GIJ,IRREGULAR=IRREGULAR,X=X,Y=Y,Z=Z,NN=NN,INTX=INTX,INTY=INTY,DERX=DERX,DERY=DERY

  ;build the polynomial basis and
  ;the basis spatial covariance matrix
  ;over the image mask pixels

  if not keyword_set(BASIS) then begin
    jnm=lonarr(3,nmodes)
    if not keyword_set(IRREGULAR) then begin
      sz=size(IMAGE)
      xyrt=COOGRID(sz[1],sz[2],scale=1)
      basis=dblarr(sz[1],sz[2],nmodes)
      j_mode=-1
      for n=0,NN do for k=0,n do begin
        j_mode=j_mode+1
        basis[*,*,j_mode]=xyrt.x^(n-k)*xyrt.y^k
        jnm[*,j_mode]=[j_mode,n-k,k]
      endfor
    endif else begin
      basis=dblarr(n_elements(X),nmodes)
      j_mode=-1
      for n=0,NN do for k=0,n do begin
        j_mode=j_mode+1
        basis[*,j_mode]=X^(n-k)*Y^k
        jnm[*,j_mode]=[j_mode,n-k,k]
      endfor
    endelse
    Gij=dblarr(nmodes,nmodes)
    if not keyword_set(IRREGULAR) then begin
      for i=0,nmodes-1 do for j=0,i do begin
        Gij[i,j]=total(basis[*,*,i]*basis[*,*,j]*MASK)
        Gij[j,i]=Gij[i,j]
      endfor
    endif else begin
      for i=0,nmodes-1 do for j=0,i do begin
        Gij[i,j]=total(basis[*,i]*basis[*,j])
        Gij[j,i]=Gij[i,j]
      endfor
    endelse
  endif
  print, '------------------'
  print, Gij[*,1]
  


  b=dblarr(nmodes)
  for i=0,nmodes-1 do
		b[i]=total(IMAGE*basis[*,*,i]*MASK)
endfor
  a=invert(Gij)#b

  ;build the model from the modes and the coefficients
  model=0
  
  for i=0,nmodes-1 do 
	model=model+a[i]*basis[*,*,i]
  endfor
  ;return
  
  
  
  return,{a:a,b:b,Gij:Gij,basis:basis,nmodes:nmodes,jnm:jnm,model:model}

end

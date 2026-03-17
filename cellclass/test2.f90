program test2
  character(len=100) :: fileout1, fileout2

  ! The original SDF string length is 33
  write(fileout1, '(A,3(A1,I1),1A4)') 'output/result/SDF/out_signed_dist','_',0,'_',0,'_',0,'.dat'
  
  ! The original heavi string length is 25
  write(fileout2, '(A,3(A1,I1),1A4)') 'output/result/G/out_heavi','_',0,'_',0,'_',0,'.dat'
  
  print *, "SDF  : '", trim(fileout1), "'"
  print *, "Heavi: '", trim(fileout2), "'"
end program test2

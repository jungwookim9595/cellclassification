program test
  character(len=100) :: fileout
  write(fileout, '(1A37,3(A1,I1),1A4)') 'output/result/G/out_heavi','_', 0, '_', 0, '_', 0, '.dat'
  print *, ">>", trim(fileout), "<<"
end program test

C --- linux_suppl.inc          
C
C     g77 does not have intrinsic sind, cosd....
C     They have to be declared external
C
C      real*8 sind, cosd, tand
C      real*8 asind, acosd, atand, atan2d
C      external sind, cosd, tand
C      external asind, acosd, atand, atan2d
C
C     linux supplemental fuctions
C
C sind
      real*8 function sind (x)
      implicit none
      real*8 x
      real*8 degrad
      parameter (degrad = 3.141592654 / 180.)

      sind = sin (degrad * x)
      return
      end
C
C cosd
      real*8 function cosd (x)
      implicit none
      real*8 x
      real*8 degrad
      parameter (degrad = 3.141592654 / 180.)

      cosd = cos (degrad * x)
      return
      end
C
C tand
      real*8 function tand (x)
      implicit none
      real*8 x
      real*8 degrad
      parameter (degrad = 3.141592654 / 180.)

      tand = tan (degrad * x)
      return
      end
C
C asind
      real*8 function asind (x)
      implicit none
      real*8 x
      real*8 raddeg
      parameter (raddeg =  180. / 3.141592654 )

      asind = raddeg * asin (x)
      return
      end
C
C acosd
      real*8 function acosd (x)
      implicit none
      real*8 x
      real*8 raddeg
      parameter (raddeg =  180. / 3.141592654 )

      acosd = raddeg * acos (x)
      return
      end
C
C atand
      real*8 function atand(x)
      implicit none
      real*8 x
      real*8 raddeg
      parameter (raddeg =  180. / 3.141592654 )

      atand = raddeg * atan (x)
      return
      end
C
C atan2d
      real*8 function atan2d (x, y)
      implicit none
      real*8 x, y
      real*8 raddeg
      parameter (raddeg =  180. / 3.141592654 )

      atan2d = raddeg * atan2 (x, y)
      return
      end
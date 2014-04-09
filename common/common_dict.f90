!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: common_dict.f90                                                   #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
!#                                                                           #
!# This program is free software; you can redistribute it and/or modify      #
!# it under the terms of the GNU General Public License as published by      #
!# the Free Software Foundation; either version 2 of the License, or (at     #
!# your option) any later version.                                           #
!#                                                                           #
!# This program is distributed in the hope that it will be useful, but       #
!# WITHOUT ANY WARRANTY; without even the implied warranty of                #
!# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
!# NON INFRINGEMENT.  See the GNU General Public License for more            #
!# details.                                                                  #
!#                                                                           #
!# You should have received a copy of the GNU General Public License         #
!# along with this program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################

!----------------------------------------------------------------------------!
!> a simple dictionary implementation for all config values
!----------------------------------------------------------------------------!
MODULE common_dict
  USE dict_common, InitDict_common => InitDict
#ifdef HAVE_NETCDF
  USE netcdf
#endif
  !--------------------------------------------------------------------------!
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE SetAttr
    MODULE PROCEDURE SetAttr1, SetAttr2, SetAttr3, SetAttr4, SetAttr5, &
        SetAttr6, SetAttr7, SetAttr7b, SetAttr8, SetAttr8b, SetAttr9, &
        SetAttr9b, SetAttr10
  END INTERFACE SetAttr
  INTERFACE GetAttr
    MODULE PROCEDURE GetAttr1, GetAttr2, GetAttr3, GetAttr4, GetAttr5, &
        GetAttr6, GetAttr7, GetAttr8, GetAttr9, GetAttr10
  END INTERFACE GetAttr
  INTERFACE GetDataType
    MODULE PROCEDURE GetDataType1, GetDataType2
  END INTERFACE GetDataType
  INTERFACE RequireKey
    MODULE PROCEDURE RequireKey0, RequireKey1, RequireKey2, RequireKey3, &
        RequireKey4, RequireKey5, RequireKey7, RequireKey8, RequireKey9, &
        RequireKey10
  END INTERFACE RequireKey
  INTERFACE CheckKey
    MODULE PROCEDURE CheckKey1, CheckKey2, CheckKey3
  END INTERFACE
!  INTERFACE ASSIGNMENT (=)
!    MODULE PROCEDURE CreateDict
!  END INTERFACE
  INTERFACE OPERATOR (/)
    MODULE PROCEDURE Assign1, Assign2, Assign3, Assign4, Assign5, Assign6, &
        Assign7, Assign8, Assign9, Assign10
  END INTERFACE
  INTERFACE AddField
    MODULE PROCEDURE AddField_2, AddField_7, AddField_8, AddField_9
  END INTERFACE
  !> \endcond
  ! constants
  INTEGER, PARAMETER   :: MAX_CHAR_LEN = 128
  INTEGER, PARAMETER   :: DICT_NONE = 0
  INTEGER, PARAMETER   :: DICT_INT  = 1
  INTEGER, PARAMETER   :: DICT_REAL = 2
  INTEGER, PARAMETER   :: DICT_CHAR = 3
  INTEGER, PARAMETER   :: DICT_BOOL = 4
  INTEGER, PARAMETER   :: DICT_REAL_ONED = 5
  INTEGER, PARAMETER   :: DICT_DIR  = 6
  INTEGER, PARAMETER   :: DICT_REAL_TWOD = 7
  INTEGER, PARAMETER   :: DICT_REAL_THREED = 8
  INTEGER, PARAMETER   :: DICT_REAL_FOURD = 9
  INTEGER, PARAMETER   :: DICT_INT_ONED = 10
  TYPE(Dict2_TYP),SAVE :: this
#define TYPE_DICT_KEY CHARACTER(LEN=MAX_CHAR_LEN)
#define TYPE_DICT_INT INTEGER
#define TYPE_DICT_REAL REAL
#define TYPE_DICT_CHAR CHARACTER(LEN=MAX_CHAR_LEN)
#define TYPE_DICT_BOOL LOGICAL
  ! common data structure
  TYPE Dict_TYP
     PRIVATE
     CHARACTER(LEN=MAX_CHAR_LEN)    :: key = ""
     INTEGER                        :: type = 0
     TYPE_DICT_INT,POINTER          :: val1 => null()
     TYPE_DICT_REAL,POINTER         :: val2 => null()
     TYPE_DICT_CHAR,POINTER         :: val3 => null()
     TYPE_DICT_BOOL,POINTER         :: val4 => null()
     TYPE_DICT_REAL,DIMENSION(:),POINTER &
                                    :: val5 => null()
     TYPE(Dict_TYP),POINTER         :: val6 => null()
     TYPE_DICT_REAL,DIMENSION(:,:),POINTER &
                                    :: val7 => null()
     TYPE_DICT_REAL,DIMENSION(:,:,:),POINTER &
                                    :: val8 => null()
     TYPE_DICT_REAL,DIMENSION(:,:,:,:),POINTER &
                                    :: val9 => null()
     TYPE_DICT_INT,DIMENSION(:),POINTER &
                                    :: val10 => null()
     TYPE(Dict_TYP),POINTER         :: next => null()
  END TYPE Dict_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Dict_TYP, &
       ! constants
       MAX_CHAR_LEN, DICT_DIR, DICT_INT, DICT_REAL, DICT_CHAR, DICT_BOOL, &
       DICT_REAL_ONED, DICT_REAL_TWOD, DICT_REAL_THREED, DICT_REAL_FOURD, &
       ! methods
       SetAttr, &
       GetAttr, &
       GetDataType, &
       GetNext, &
       GetKey, &
       GetVal1, GetVal2, GetVal3, GetVal4, GetVal5, GetVal6, GetVal7, &
       GetVal8, GetVal9, GetVal10, &
       OPERATOR(/), &
!       OPERATOR(+), &
!       ASSIGNMENT(=), &
       AddField, &
       Dict, &
       RequireKey, &
       HasKey, &
!       NumSubKeys, &
       CheckKey, &
       CopyHierarchy, &
#ifdef HAVE_NETCDF
       SaveDict, &       
       LoadDict, &
#endif
       PrintNode, &
       PrintDict, &
       InitDict, &
       DeleteDict
  !--------------------------------------------------------------------------!

CONTAINS
  SUBROUTINE InitDict()
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, PARAMETER :: type = 0
    CHARACTER(LEN=32), PARAMETER :: name = "Dictionary"
    !------------------------------------------------------------------------!
    CALL InitDict_common(this, type, name)
  END SUBROUTINE InitDict

  SUBROUTINE AppendDict(root, appendix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: root
    TYPE(Dict_TYP),TARGET :: appendix
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: node,app,tmp
    TYPE_DICT_CHAR  :: k
    INTEGER :: i

    !------------------------------------------------------------------------!
    INTENT(IN)      :: appendix
    !------------------------------------------------------------------------!
    IF(ASSOCIATED(root)) THEN
       app => appendix
       k = GetKey(app)
       node => root
       DO 
         IF (GetKey(node) .eq. k) THEN
           !print *, "beep ", TRIM(GetKey(node)),TRIM(k)
           node => node%val6
           app => app%val6
           k = TRIM(GetKey(app))
         ELSE IF (.NOT.ASSOCIATED(node%next)) THEN
           EXIT            
         ELSE
           node => node%next
         END IF
       END DO
       node%next => app
    ELSE
        root => appendix
    END IF
  END SUBROUTINE AppendDict

!  SUBROUTINE CreateDict(root, appendix)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Dict_TYP),POINTER &
!                    :: root, appendix
!    !------------------------------------------------------------------------!
!    root => appendix
!  END SUBROUTINE CreateDict
  FUNCTION  Dict(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,&
                 n18,n19,n20) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: res,n1
    TYPE(Dict_TYP),POINTER,OPTIONAL :: n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,&
                                       n13,n14,n15,n16,n17,n18,n19,n20
    !------------------------------------------------------------------------!
    res => n1
    IF(PRESENT(n2)) CALL AppendDict(res, n2)
    IF(PRESENT(n3)) CALL AppendDict(res, n3)
    IF(PRESENT(n4)) CALL AppendDict(res, n4)
    IF(PRESENT(n5)) CALL AppendDict(res, n5)
    IF(PRESENT(n6)) CALL AppendDict(res, n6)
    IF(PRESENT(n7)) CALL AppendDict(res, n7)
    IF(PRESENT(n8)) CALL AppendDict(res, n8)
    IF(PRESENT(n9)) CALL AppendDict(res, n9)
    IF(PRESENT(n10)) CALL AppendDict(res, n10)
    IF(PRESENT(n11)) CALL AppendDict(res, n11)
    IF(PRESENT(n12)) CALL AppendDict(res, n12)
    IF(PRESENT(n13)) CALL AppendDict(res, n13)
    IF(PRESENT(n14)) CALL AppendDict(res, n14)
    IF(PRESENT(n15)) CALL AppendDict(res, n15)
    IF(PRESENT(n16)) CALL AppendDict(res, n16)
    IF(PRESENT(n17)) CALL AppendDict(res, n17)
    IF(PRESENT(n18)) CALL AppendDict(res, n18)
    IF(PRESENT(n19)) CALL AppendDict(res, n19)
    IF(PRESENT(n20)) CALL AppendDict(res, n20)
  END FUNCTION Dict


  FUNCTION Assign1(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_INT   :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION

  FUNCTION Assign2(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL  :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign2

  FUNCTION Assign3(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    CHARACTER(LEN=*):: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign3

  FUNCTION Assign4(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_BOOL  :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign4

  FUNCTION Assign5(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign5

  FUNCTION Assign6(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE(Dict_TYP),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign6

  FUNCTION Assign7(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign7

  FUNCTION Assign8(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:,:),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign8

  FUNCTION Assign9(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:,:,:),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign9

  FUNCTION Assign10(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    TYPE_DICT_INT,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign10

  FUNCTION IsType(root, type) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    INTEGER         :: type
    LOGICAL         :: res
    !------------------------------------------------------------------------!
    INTENT(IN)      :: type
    !------------------------------------------------------------------------!
    res = (root%type.EQ.type)
  END FUNCTION IsType


  RECURSIVE FUNCTION TouchKey(root, key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res, dir
    CHARACTER(LEN=*):: key
    TYPE_DICT_CHAR  :: k
    INTEGER         :: i
    !------------------------------------------------------------------------!
    INTENT(IN)         :: key
    !------------------------------------------------------------------------!
    k = key
    IF(k(1:1).EQ.'/') k = k(2:)
    i = SCAN(k, '/', .TRUE.)
    !print *,"key=",k
    dir => FindDir(root, k)
    !print *,"Dir?: ", ASSOCIATED(dir)
    IF((i.NE.0).AND.(.NOT.ASSOCIATED(dir))) THEN
        dir => TouchKey(root, k(1:i-1))
        ALLOCATE(res)
        res%key = k(i+1:)
        CALL SetAttr(dir, dir%key, res)
    ELSE
        res => FindKey(dir, k)
        ! Didn't find key? create new node!
        IF(.NOT.ASSOCIATED(res)) THEN
            ALLOCATE(res)
            res%key = k(i+1:)
            CALL AppendDict(dir, res)
        END IF
    END IF
    IF(.NOT.ASSOCIATED(root)) root => dir
  END FUNCTION TouchKey

  RECURSIVE FUNCTION SetAttr0(root, key, type, sh) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    INTEGER         :: type
    INTEGER,DIMENSION(:),OPTIONAL&
                    :: sh
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, type
    !------------------------------------------------------------------------!
    res => TouchKey(root, key)
    CALL TouchType(res, type, sh)
  END FUNCTION SetAttr0

  SUBROUTINE SetAttr1(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_INT   :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_INT)
    res%val1 = val
  END SUBROUTINE SetAttr1

  SUBROUTINE SetAttr2(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL  :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL)
    res%val2 = val
  END SUBROUTINE SetAttr2

  SUBROUTINE SetAttr3(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    CHARACTER(LEN=*):: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_CHAR)
    res%val3 = val
  END SUBROUTINE SetAttr3

  SUBROUTINE SetAttr4(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_BOOL  :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    CALL Error(this, "SetAttr4", "LOGICAL not supported anymore!")
    res => SetAttr0(root,key,DICT_BOOL)
    res%val4 = val
  END SUBROUTINE SetAttr4

  SUBROUTINE SetAttr5(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_ONED,SHAPE(val))
    res%val5 = val
  END SUBROUTINE SetAttr5

  RECURSIVE SUBROUTINE SetAttr6(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE(Dict_TYP),TARGET &
                    :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_DIR)
    res%val6 => val
  END SUBROUTINE SetAttr6

  SUBROUTINE SetAttr7(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_TWOD,SHAPE(val))
!    res%val7 = val
    res%val7 => val
  END SUBROUTINE SetAttr7

  SUBROUTINE SetAttr7b(root, key, val,bn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    INTEGER,DIMENSION(2),INTENT(IN) :: bn
    TYPE_DICT_REAL,DIMENSION(bn(1):,bn(2):),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_TWOD,SHAPE(val))
!    res%val7 = val
    res%val7 => val
  END SUBROUTINE SetAttr7b

  SUBROUTINE SetAttr8(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:,:),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_THREED,SHAPE(val))
!    res%val8 = val
    res%val8 => val
  END SUBROUTINE SetAttr8

  SUBROUTINE SetAttr8b(root, key, val,bn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    INTEGER,DIMENSION(3),INTENT(IN) :: bn
    TYPE_DICT_REAL,DIMENSION(bn(1):,bn(2):,bn(3):),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_THREED,SHAPE(val))
!    res%val8 = val
    res%val8 => val
  END SUBROUTINE SetAttr8b

  SUBROUTINE SetAttr9(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:,:,:),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_FOURD,SHAPE(val))
!    res%val9 = val
    res%val9 => val
  END SUBROUTINE SetAttr9

  SUBROUTINE SetAttr9b(root, key, val, bn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    INTEGER,DIMENSION(4),INTENT(IN) :: bn
    TYPE_DICT_REAL,DIMENSION(bn(1):,bn(2):,bn(3):,bn(4):),TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_REAL_FOURD,SHAPE(val))
!    res%val9 = val
    res%val9 => val
  END SUBROUTINE SetAttr9b

  SUBROUTINE SetAttr10(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    CHARACTER(LEN=*):: key
    TYPE_DICT_INT,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    res => SetAttr0(root,key,DICT_INT_ONED,SHAPE(val))
    res%val10 = val
  END SUBROUTINE SetAttr10

  RECURSIVE SUBROUTINE TouchType(node, type, sh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: node
    INTEGER         :: type
    INTEGER, DIMENSION(:), OPTIONAL &
                    :: sh
    !------------------------------------------------------------------------!
    INTENT(IN)      :: type, sh
    !------------------------------------------------------------------------!
    
    IF(node%type.NE.type) THEN
        SELECT CASE(node%type)
        CASE(1)
            DEALLOCATE(node%val1)
            NULLIFY(node%val1)
        CASE(2)
            DEALLOCATE(node%val2)
            NULLIFY(node%val2)
        CASE(3)
            DEALLOCATE(node%val3)
            NULLIFY(node%val3)
        CASE(4)
            DEALLOCATE(node%val4)
            NULLIFY(node%val4)
        CASE(5)
            DEALLOCATE(node%val5)
            NULLIFY(node%val5)
        CASE(6)
            CALL DeleteDict(node%val6)
            NULLIFY(node%val6)
        CASE(7)
!            DEALLOCATE(node%val7)
            NULLIFY(node%val7)
        CASE(8)
!            DEALLOCATE(node%val8)
            NULLIFY(node%val8)
        CASE(9)
!            DEALLOCATE(node%val9)
            NULLIFY(node%val9)
        CASE(10)
            DEALLOCATE(node%val10)
            NULLIFY(node%val10)
        END SELECT
        node%type = type
        SELECT CASE(node%type)
        CASE(1)
            ALLOCATE(node%val1)
        CASE(2)
            ALLOCATE(node%val2)
        CASE(3)
            ALLOCATE(node%val3)
        CASE(4)
            ALLOCATE(node%val4)
        CASE(5)
            ALLOCATE(node%val5(sh(1)))
        CASE(6)
            ! Do nothing
        CASE(7)
!            ALLOCATE(node%val7(sh(1),sh(2)))
        CASE(8)
!            ALLOCATE(node%val8(sh(1),sh(2),sh(3)))
        CASE(9)
!            ALLOCATE(node%val9(sh(1),sh(2),sh(3),sh(4)))
        CASE(10)
            ALLOCATE(node%val10(sh(1)))
        END SELECT
    END IF

  END SUBROUTINE TouchType

RECURSIVE SUBROUTINE PrintNode(node, prefix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: node
    CHARACTER(LEN=10*MAX_CHAR_LEN)  :: val
    CHARACTER(LEN=*),OPTIONAL &
                    :: prefix
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!

    !print *,"type: ", node%type
    SELECT CASE(node%type)
    CASE(DICT_NONE)
        val = ''
    CASE(DICT_INT)
        WRITE(val,'(I8)')node%val1
    CASE(DICT_REAL)
        WRITE(val,*)node%val2
    CASE(DICT_CHAR)
        WRITE(val,'(A)')node%val3
    CASE(DICT_BOOL)
        WRITE(val,*)node%val4
    CASE(DICT_REAL_ONED)
!        WRITE(val,*)node%val5
        WRITE(val,*)"1D real array, shape: ",SHAPE(node%val5)
    CASE(DICT_DIR)
        IF(PRESENT(prefix)) THEN
            val = TRIM(prefix)//'/'//TRIM(node%key)
        ELSE
            val = '/'//TRIM(node%key)
        END IF
        CALL PrintDict(node%val6,val)
    CASE(DICT_REAL_TWOD)
        !WRITE(val,*)node%val7
        WRITE(val,*)"2D real array, shape: ",SHAPE(node%val7)
    CASE(DICT_REAL_THREED)
        !WRITE(val,*)node%val8
        WRITE(val,*)"3D real array, shape: ",SHAPE(node%val8)
    CASE(DICT_REAL_FOURD)
        !WRITE(val,*)node%val9
        WRITE(val,*)"4D real array, shape: ",SHAPE(node%val9)
    CASE(DICT_INT_ONED)
        WRITE(val,*)"1D int array, shape: ",SHAPE(node%val10)
    END SELECT
    IF(node%type.NE.6) THEN
        IF(PRESENT(prefix)) THEN
            CALL Info(this,TRIM(prefix)//'/'//TRIM(node%key)//': '//TRIM(val))
        ELSE
            CALL Info(this,TRIM(node%key)//' '//TRIM(val))
        END IF
    END IF
  END SUBROUTINE PrintNode

  RECURSIVE SUBROUTINE PrintDict(root, prefix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node
    CHARACTER(LEN=*),OPTIONAL &
                    :: prefix
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    node => root
    DO WHILE(ASSOCIATED(node))
        CALL PrintNode(node, prefix)
        node => node%next
    END DO
  END SUBROUTINE PrintDict

  FUNCTION FindKey(root, key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res, node
    CHARACTER(LEN=*):: key
    CHARACTER(LEN=MAX_CHAR_LEN) :: k
    INTEGER         :: i
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    NULLIFY(res)
    node => root
    i = SCAN(key, '/',.TRUE.)
    k = key
    k = k(i+1:)
    DO WHILE(ASSOCIATED(node).AND.(.NOT.ASSOCIATED(res)))
        !print *,"Findkey1=",node%key
        !print *,"Findkey2=",k
        IF(node%key.EQ.k) THEN
            res => node
        ELSE
            node => node%next
        END IF
    END DO
        
  END FUNCTION FindKey

!  RECURSIVE FUNCTION NumSubKeys(root, key) RESULT(res)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Dict_TYP),POINTER &
!                    :: root, node
!    CHARACTER(LEN=*):: key
!    CHARACTER(LEN=MAX_CHAR_LEN) :: k
!    INTEGER         :: res,i
!    !------------------------------------------------------------------------!
!    INTENT(IN)      :: key
!    !------------------------------------------------------------------------!
!    node => root
!    i = SCAN(key, '/',.TRUE.)
!    k = key
!    k = k(i+1:)
!    res = 0
!    DO WHILE(ASSOCIATED(node))
!        IF(node%key.EQ.k) res = res +1
!        IF (IsType(node,DICT_DIR)) res = res + NumSubKeys(node%val6, k)
!        node => node%next
!    END DO
!  END FUNCTION NumSubKeys

  FUNCTION HasKey(root, key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=*):: key
    LOGICAL         :: res
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    !print *,"<<<<<<<<<<<<",TRIM(key)
    !CALL PrintDict(root)
    !print *,"------------"
    dir => FindDir(root,key)
    !CALL PrintDict(dir)
    !print *,">>>>>>>>>>>"
    res = ASSOCIATED(FindKey(dir,key))
  END FUNCTION HasKey

  FUNCTION CheckKey1(root, key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=*):: key
    LOGICAL         :: res
    INTEGER         :: val, val2
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res = HasKey(root, key)
    IF(res) THEN
        res = .FALSE.
        CALL GetAttr(root, key, val2)
        IF(val.EQ.val2) THEN
            res = .TRUE.
        END IF
    END IF    
  END FUNCTION CheckKey1

  FUNCTION CheckKey2(root, key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=*):: key
    LOGICAL         :: res
    REAL            :: val, val2
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res = HasKey(root, key)
    IF(res) THEN
        res = .FALSE.
        CALL GetAttr(root, key, val2)
        IF(val.EQ.val2) THEN
            res = .TRUE.
        END IF
    END IF    
  END FUNCTION CheckKey2

  FUNCTION CheckKey3(root, key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=*):: key
    LOGICAL         :: res
    CHARACTER(LEN=MAX_CHAR_LEN) :: val, val2
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    res = HasKey(root, key)
    IF(res) THEN
        res = .FALSE.
        CALL GetAttr(root, key, val2)
        IF(val.EQ.val2) THEN
            res = .TRUE.
        END IF
    END IF    
  END FUNCTION CheckKey3

  FUNCTION GetKey(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=MAX_CHAR_LEN):: res
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    res = root%key
  END FUNCTION GetKey

  RECURSIVE SUBROUTINE CopyHierarchy(root, outdir)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, outdir, dir,temp,subdir
    CHARACTER(LEN=MAX_CHAR_LEN) :: key
    !------------------------------------------------------------------------!
    dir => root
    DO WHILE(ASSOCIATED(dir))
       IF (IsType(dir,DICT_DIR)) THEN
          key = TRIM(GetKey(dir))
          temp => Dict("." / 0)
          CALL SetAttr(outdir, key, temp)
          CALL GetAttr(dir, key, subdir)
          CALL CopyHierarchy(subdir,temp)
       END IF
       dir => GetNext(dir)
    END DO
  END SUBROUTINE CopyHierarchy

  RECURSIVE FUNCTION FindDir(root, key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res, node
    CHARACTER(LEN=*):: key
    TYPE_DICT_CHAR  :: k
    INTEGER         :: i
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------! 
    k = key
    IF(k(1:1).EQ.'/') k = k(2:)
    i = SCAN(k, '/')
    IF((i.NE.0).AND.ASSOCIATED(root)) THEN
        !print *,"Search: ", k(1:i-1)
        CALL GetAttr(root, k(1:i-1), node)
        !print *,"Found: ", ASSOCIATED(node)
        res => FindDir(node, k(i+1:))
    ELSE
        res => root
    END IF
  END FUNCTION FindDir

  RECURSIVE SUBROUTINE GetAttr0(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, dir, res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      ::  key
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    res => FindKey(dir, key)
  END SUBROUTINE GetAttr0

  SUBROUTINE GetAttr1(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_INT   :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    INTENT(INOUT)   :: res
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_INT)) &
        CALL Error(this, "GetAttr1", "Key '"//TRIM(key)//"' has wrong type")
        res = node%val1
    ELSE
        CALL Error(this, "GetAttr1", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr1

  SUBROUTINE GetAttr2(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_REAL  :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      ::  key
    INTENT(INOUT)   :: res
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_REAL)) &
        CALL Error(this, "GetAttr2", "Key '"//TRIM(key)//"' has wrong type")
        res = node%val2
    ELSE
        CALL Error(this, "GetAttr2", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr2
  
  SUBROUTINE GetAttr3(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_CHAR  :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    INTENT(INOUT)   :: res
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_CHAR)) &
        CALL Error(this, "GetAttr3", "Key '"//TRIM(key)//"' has wrong type")
        res = node%val3
    ELSE
        CALL Error(this, "GetAttr3", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr3
  
  SUBROUTINE GetAttr4(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_BOOL  :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    INTENT(INOUT)   :: res
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_BOOL)) &
        CALL Error(this, "GetAttr4", "Key '"//TRIM(key)//"' has wrong type")
        res = node%val4
    ELSE
        CALL Error(this, "GetAttr4", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr4

  SUBROUTINE GetAttr5(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_REAL, DIMENSION(:) :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    INTENT(INOUT)   :: res
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_REAL_ONED)) &
        CALL Error(this, "GetAttr5", "Key '"//TRIM(key)//"' has wrong type")
      IF (SIZE(res) .NE. SIZE(node%val5)) &
        CALL Error(this, "GetAttr5", "Key '"//TRIM(key)//"' has wrong dimension")
        res = node%val5
    ELSE
        CALL Error(this, "GetAttr5", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr5

  SUBROUTINE GetAttr6(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE(Dict_TYP),POINTER &
                    :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_DIR)) &
        CALL Error(this, "GetAttr6", "Key '"//TRIM(key)//"' has wrong type")
        res => node%val6
    END IF
  END SUBROUTINE GetAttr6

  SUBROUTINE GetAttr7(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_REAL, DIMENSION(:,:),POINTER :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_REAL_TWOD)) &
        CALL Error(this, "GetAttr7", "Key '"//TRIM(key)//"' has wrong type")
        res => node%val7
    ELSE
        CALL Error(this, "GetAttr7", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr7

  SUBROUTINE GetAttr8(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_REAL, DIMENSION(:,:,:),POINTER :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_REAL_THREED)) &
        CALL Error(this, "GetAttr8", "Key '"//TRIM(key)//"' has wrong type")
        res => node%val8
    ELSE
        CALL Error(this, "GetAttr8", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr8

  SUBROUTINE GetAttr9(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_REAL, DIMENSION(:,:,:,:),POINTER :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_REAL_FOURD)) &
        CALL Error(this, "GetAttr9", "Key '"//TRIM(key)//"' has wrong type")
        res => node%val9
    ELSE
        CALL Error(this, "GetAttr9", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr9

  SUBROUTINE GetAttr10(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node => null()
    TYPE_DICT_INT, DIMENSION(:) :: res
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    INTENT(INOUT)   :: res
    !------------------------------------------------------------------------!
    CALL GetAttr0(root, key, node)
    IF(ASSOCIATED(node)) THEN
      IF(.NOT.isType(node, DICT_INT_ONED)) &
        CALL Error(this, "GetAttr10", "Key '"//TRIM(key)//"' has wrong type")
      IF (SIZE(res) .NE. SIZE(node%val10)) &
        CALL Error(this, "GetAttr10", "Key '"//TRIM(key)//"' has wrong dimension")
        res = node%val10
    ELSE
        CALL Error(this, "GetAttr10", "Couldn't find key '"//TRIM(key)//"'.")
    END IF
  END SUBROUTINE GetAttr10

  SUBROUTINE RequireKey0(root, key)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=*):: key
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    IF(.NOT.HasKey(root,key)) THEN
        CALL Error(this, "RequireKey0","A required argument with key='"//&
                                        TRIM(key)//"' is undefined.")
    END IF
  END SUBROUTINE RequireKey0

  SUBROUTINE RequireKey1(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_INT   :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_INT)) THEN
        CALL Warning(this, "RequireKey1", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey1

  SUBROUTINE RequireKey2(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL  :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_REAL)) THEN
        CALL Warning(this, "RequireKey2", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey2

  SUBROUTINE RequireKey3(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_CHAR  :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_CHAR)) THEN
        CALL Warning(this, "RequireKey3", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey3

  SUBROUTINE RequireKey4(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_BOOL  :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_BOOL)) THEN
        CALL Warning(this, "RequireKey4", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey4

  SUBROUTINE RequireKey5(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key, val
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_REAL_ONED)) THEN
        CALL Warning(this, "RequireKey5", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey5

  SUBROUTINE RequireKey7(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:),POINTER :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_REAL_TWOD)) THEN
        CALL Warning(this, "RequireKey7", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey7
  
  SUBROUTINE RequireKey8(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:,:),POINTER :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_REAL_THREED)) THEN
        CALL Warning(this, "RequireKey8", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey8

  SUBROUTINE RequireKey9(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_REAL,DIMENSION(:,:,:,:),POINTER :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_REAL_FOURD)) THEN
        CALL Warning(this, "RequireKey9", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey9

  SUBROUTINE RequireKey10(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, tmp
    CHARACTER(LEN=*):: key
    TYPE_DICT_INT,DIMENSION(:) :: val
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)      :: key
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    tmp => FindKey(dir, key)
    IF(.NOT.ASSOCIATED(tmp)) THEN
        CALL SetAttr(root, key, val)
    ELSE IF(.NOT.isType(tmp, DICT_INT_ONED)) THEN
        CALL Warning(this, "RequireKey10", "Key '"//TRIM(tmp%key)//"' has wrong type")
        CALL SetAttr(root, key, val)
    END IF
  END SUBROUTINE RequireKey10

  RECURSIVE SUBROUTINE DeleteNode(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    INTEGER         :: status
    CHARACTER(LEN=32) :: str = ""
    !------------------------------------------------------------------------!
    CALL TouchType(root, 0)
    DEALLOCATE(root, stat=status)
    ! nesh: Error code 195 means deallocation of a unassociated pointer.
    ! Checking with ASSOCIATED says it is indeed associated. So ignore error
    ! 195
    IF((status.NE.0).AND.(status.NE.195)) THEN
       WRITE(str, *)status
       CALL Error(this,"DeleteNode", "status = "//TRIM(str))
    END IF
    NULLIFY(root)
  END SUBROUTINE DeleteNode

  RECURSIVE SUBROUTINE DeleteDict(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node, next
    !------------------------------------------------------------------------!
    node => root
    DO WHILE(ASSOCIATED(node))
        next => node%next
        IF (IsType(node,DICT_DIR)) THEN
           CALL DeleteDict(node%val6)
        ELSE
           CALL DeleteNode(node)
        END IF
        node => next
    END DO
  END SUBROUTINE DeleteDict

  FUNCTION GetNext(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, res
    !------------------------------------------------------------------------!
    res => root%next
  END FUNCTION GetNext

  FUNCTION GetDataType1(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    INTEGER         :: res
    !------------------------------------------------------------------------!
    res = root%type
  END FUNCTION GetDataType1

  FUNCTION GetDataType2(root,key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root,dir
    CHARACTER(LEN=*):: key
    INTEGER         :: res
    !------------------------------------------------------------------------!
    dir => FindDir(root, key)
    dir => FindKey(dir, key)
    res = dir%type
  END FUNCTION GetDataType2

  FUNCTION GetVal1(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    INTEGER,POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_INT)) &
        CALL Error(this, "GetVal1", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val1
  END FUNCTION GetVal1

  FUNCTION GetVal2(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    REAL,POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_REAL)) &
        CALL Error(this, "GetVal2", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val2
  END FUNCTION GetVal2

  FUNCTION GetVal3(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    CHARACTER(LEN=MAX_CHAR_LEN),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_CHAR)) &
        CALL Error(this, "GetVal3", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val3
  END FUNCTION GetVal3

  FUNCTION GetVal4(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    LOGICAL,POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_BOOL)) &
        CALL Error(this, "GetVal4", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val4
  END FUNCTION GetVal4

  FUNCTION GetVal5(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    REAL,DIMENSION(:),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_REAL_ONED)) &
        CALL Error(this, "GetVal5", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val5
  END FUNCTION GetVal5

  FUNCTION GetVal6(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    TYPE(Dict_TYP),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_DIR)) &
        CALL Error(this, "GetVal6", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val6
  END FUNCTION GetVal6

  FUNCTION GetVal7(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    REAL,DIMENSION(:,:),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_REAL_TWOD)) &
        CALL Error(this, "GetVal7", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val7
  END FUNCTION GetVal7

  FUNCTION GetVal8(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    REAL,DIMENSION(:,:,:),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_REAL_THREED)) &
        CALL Error(this, "GetVal8", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val8
  END FUNCTION GetVal8

  FUNCTION GetVal9(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    REAL,DIMENSION(:,:,:,:),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_REAL_FOURD)) &
        CALL Error(this, "GetVal9", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val9
  END FUNCTION GetVal9

  FUNCTION GetVal10(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root
    TYPE_DICT_INT,DIMENSION(:),POINTER :: res
    !------------------------------------------------------------------------!
    IF(.NOT.isType(root, DICT_INT_ONED)) &
        CALL Error(this, "GetVal10", "Key '"//TRIM(root%key)//"' has wrong type")
    res => root%val10
  END FUNCTION GetVal10

  SUBROUTINE AddField_2(IO, key, value, attr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: IO
    CHARACTER(LEN=*)            :: key
    REAL,POINTER                :: value
    TYPE(Dict_TYP),POINTER,OPTIONAL :: attr
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: dic => null()
    !------------------------------------------------------------------------!
    INTENT(IN)  :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(attr)) THEN
        dic => attr
    END IF
    CALL SetAttr(dic, "value", value)
    CALL SetAttr(IO, key, dic)
  END SUBROUTINE AddField_2

  SUBROUTINE AddField_7(IO, key, field, attr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: IO
    CHARACTER(LEN=*)            :: key
    REAL,DIMENSION(:,:),POINTER :: field
    TYPE(Dict_TYP),POINTER,OPTIONAL :: attr
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: dic => null()
    !------------------------------------------------------------------------!
    INTENT(IN)  :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(attr)) THEN
        dic => attr
    END IF
    CALL SetAttr(dic, "value", field,LBOUND(field))
    CALL SetAttr(IO, key, dic)
  END SUBROUTINE AddField_7

  SUBROUTINE AddField_8(IO, key, field, attr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: IO
    CHARACTER(LEN=*)            :: key
    REAL,DIMENSION(:,:,:),POINTER :: field
    TYPE(Dict_TYP),POINTER,OPTIONAL :: attr
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: dic => null()
    !------------------------------------------------------------------------!
    INTENT(IN)  :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(attr)) THEN
        dic => attr
    END IF
    CALL SetAttr(dic, "value", field,LBOUND(field))
    CALL SetAttr(IO, key, dic)
  END SUBROUTINE AddField_8

  SUBROUTINE AddField_9(IO, key, field, attr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: IO
    CHARACTER(LEN=*)            :: key
    REAL,DIMENSION(:,:,:,:),POINTER :: field
    TYPE(Dict_TYP),POINTER,OPTIONAL :: attr
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER      :: dic => null()
    !------------------------------------------------------------------------!
    INTENT(IN)  :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(attr)) THEN
        dic => attr
    END IF
    CALL SetAttr(dic, "value", field,LBOUND(field))
    CALL SetAttr(IO, key, dic)
  END SUBROUTINE AddField_9

#ifdef HAVE_NETCDF
  RECURSIVE FUNCTION SaveNode(node, rootid, ncid) RESULT(status)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: node
    INTEGER         :: ncid, rootid
    INTEGER         :: status
    !------------------------------------------------------------------------!
    INTEGER         :: grpid, varid, jnum, inum
    INTEGER         :: twoddims(2)
    !------------------------------------------------------------------------!
    INTENT(IN)      :: ncid, rootid
    !------------------------------------------------------------------------!
    status = NF90_NOERR
    SELECT CASE(node%type)
    CASE(DICT_NONE)
    CASE(DICT_INT)
        status = NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(node%key),node%val1)
    CASE(DICT_REAL)
        status = NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(node%key),node%val2)
    CASE(DICT_CHAR)
        status = NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(node%key),node%val3)
    CASE(DICT_BOOL)
        CALL Warning(this, "SaveNode", "Saving of logical types is not "&
                //"supported")
        !status = NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(node%key),node%val4)
    CASE(DICT_REAL_ONED)
        status = NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(node%key),node%val5)
    CASE(DICT_DIR)
        status = NF90_DEF_GRP(ncid, TRIM(node%key), grpid)
        IF(status.EQ.NF90_NOERR) &
            status = SaveDict(node%val6, rootid, grpid)
    CASE(DICT_REAL_TWOD)
        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQ_DIMID(ncid, "jnum", twoddims(1))
        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQ_DIMID(ncid, "inum", twoddims(2))
        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQUIRE_DIMENSION(ncid, twoddims(1), len = jnum)
        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQUIRE_DIMENSION(ncid, twoddims(2), len = inum)
        !IF(status.EQ.NF90_NOERR) &
        !    status = NF90_ENDDEF(rootid)
        !IF(status.EQ.NF90_NOERR) &
        !    status = NF90_REDEF(rootid)
        IF(status.EQ.NF90_NOERR) &
            status = NF90_DEF_VAR(ncid, TRIM(node%key), NF90_REAL8,&
                                  twoddims, varid)
        !IF(status.EQ.NF90_NOERR) &
        !    status = NF90_ENDDEF(rootid)
        IF(status.EQ.NF90_NOERR) &
            status = NF90_PUT_VAR(ncid, varid, node%val7, &
                                  map=(/ inum, 1 /), &
                                  count=(/jnum, inum/), &
                                  start=(/1, 1/))
    CASE(DICT_REAL_THREED)
        !status = NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(node%key),node%val8)
    CASE DEFAULT
        CALL Error(this, "SaveNode", "type error")
        stop
    END SELECT

  END FUNCTION SaveNode

  RECURSIVE FUNCTION SaveDict(root, rootid, ncid) RESULT(status)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER &
                    :: root, node
    INTEGER,OPTIONAL:: ncid
    INTEGER         :: rootid
    INTEGER         :: status
    !------------------------------------------------------------------------!
    INTENT(IN)      :: ncid, rootid
    !------------------------------------------------------------------------!
    status = NF90_NOERR
    node => root
    DO WHILE(ASSOCIATED(node))
        IF(status.EQ.NF90_NOERR) THEN
            IF(PRESENT(ncid)) THEN
                status = SaveNode(node, rootid, ncid)
            ELSE
                status = SaveNode(node, rootid, rootid)
            END IF
        END IF
        node => node%next
    END DO

  END FUNCTION SaveDict

  FUNCTION LoadDict(ncid,config) RESULT(status)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER             :: ncid
    TYPE(Dict_TYP), POINTER &
                        :: config
    !------------------------------------------------------------------------!
    INTEGER             :: nGrps, i, j, nVars, nAtts, nDims, xtype, n, &
                           jnum, dims(3), k
    INTEGER             :: status
    INTEGER             :: intval
    DOUBLE PRECISION    :: realval
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: realonedval
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: realtwodval
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: realthreedval
    INTEGER,DIMENSION(100)&
                        :: ncids, buf
    INTEGER,DIMENSION(NF90_MAX_VAR_DIMS)&
                        :: dimids
    CHARACTER(LEN=NF90_MAX_NAME)&
                        :: grpname,attname,varname,strval,path,key
    !------------------------------------------------------------------------!
    INTENT(IN)          :: ncid
    !------------------------------------------------------------------------!
    status = NF90_INQ_GRPS(ncid, nGrps, ncids)
    
    !Add first sublevel groups
    DO i=1,nGrps
        n = 0
        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQ_GRPS(ncids(i), n, buf)
        DO j=1,n
            ncids(nGrps+j) = buf(j)
        END DO
        nGrps = nGrps + n
    END DO

    DO i=1,nGrps 
        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQ_GRPNAME(ncids(i), grpname)

        !print *,grpname

        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQUIRE(ncids(i),nDims,nVars,nAtts)

        !print *,nAtts

        IF(status.EQ.NF90_NOERR) &
            status = NF90_INQ_GRPNAME_FULL(ncids(i), n, path)


        DO j=1,nAtts
            IF(status.EQ.NF90_NOERR) &
                status = NF90_INQ_ATTNAME(ncids(i), NF90_GLOBAL, j, attname)
            !print *,path
            !print *,attname

            IF(status.EQ.NF90_NOERR) &
                status = NF90_INQUIRE_ATTRIBUTE(ncids(i),NF90_GLOBAL,attname, &
                                                xtype, n)
            ! xtype is one of NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT,
            ! NF90_FLOAT, and NF90_DOUBLE. 

            strval=''
            key = TRIM(path)//'/'//TRIM(attname)
            IF(status.EQ.NF90_NOERR) THEN
                SELECT CASE(xtype)
                CASE(NF90_INT)
                    status = NF90_GET_ATT(ncids(i), NF90_GLOBAL, attname, intval)
                    WRITE(strval,'(I10)')intval
                    CALL SetAttr(config, key, intval)
                    !print *,intval
                CASE(NF90_INT64)
                    status = NF90_GET_ATT(ncids(i), NF90_GLOBAL, attname, intval)
                    WRITE(strval,'(I10)')intval
                    CALL SetAttr(config, key, intval)
                    !print *,intval
                CASE(NF90_DOUBLE)
                    IF(n.EQ.1) THEN
                        status = NF90_GET_ATT(ncids(i), NF90_GLOBAL, attname, realval)
                        WRITE(strval,'(E11.5)')realval
                        CALL SetAttr(config, key, realval)
                    ELSE
                        ALLOCATE(realonedval(n))
                        status = NF90_GET_ATT(ncids(i), NF90_GLOBAL, attname, &
                                              realonedval)
                        WRITE(strval,'(E11.5)')realonedval(1)
                        CALL SetAttr(config, key, realonedval)
                        DEALLOCATE(realonedval)
                    END IF
                    !print *,realval
                CASE(NF90_CHAR)
                    status = NF90_GET_ATT(ncids(i), NF90_GLOBAL, attname, strval)
                    CALL SetAttr(config, key, strval)
                    !print *,strval
                CASE DEFAULT
                    CALL Warning(this, "LoadDict_netcdf", &
                        "Unknown data type while reading an attribute.")
                END SELECT
            END IF 
            !print *,TRIM(path)//'/'//TRIM(attname)//': '//TRIM(strval)

        END DO

        DO j=1,nVars
            IF(status.EQ.NF90_NOERR) &
                status = NF90_INQUIRE_VARIABLE(ncids(i), j, varname, xtype, nDims, dimids)
            key = TRIM(path)//'/'//TRIM(varname)
            DO k=1,nDims
               IF(status.EQ.NF90_NOERR) &
                    status = NF90_INQUIRE_DIMENSION(ncids(i), dimids(k), len = dims(k))
            END DO
            SELECT CASE(nDims)
            CASE(1)
                ALLOCATE(realonedval(dims(3)))
                IF(status.EQ.NF90_NOERR) &
                    status = NF90_GET_VAR(ncids(i), j, realonedval, &
                                          map=(/ 1 /), &
                                          count=(/dims(1)/))
            CALL SetAttr(config, key, realonedval) 
            DEALLOCATE(realonedval)
            CASE(2)
                ALLOCATE(realtwodval(dims(2),dims(1)))
                IF(status.EQ.NF90_NOERR) &
                    status = NF90_GET_VAR(ncids(i), j, realtwodval, &
                                          map=(/dims(2), 1/),&
                                          count=(/dims(1),dims(2)/))
                CALL SetAttr(config, key, realtwodval) 
                DEALLOCATE(realtwodval)
            CASE(3)
                ALLOCATE(realthreedval(dims(3),dims(2),dims(1)))
                IF(status.EQ.NF90_NOERR) &
                    status = NF90_GET_VAR(ncids(i), j, realthreedval, &
                                          map=(/dims(3)*dims(2), dims(2), 1/), &
                                          count=(/dims(1),dims(2),dims(3)/))
                CALL SetAttr(config, key, realthreedval) 
                DEALLOCATE(realthreedval)
            END SELECT

        END DO
    END DO

  END FUNCTION LoadDict
#endif

END MODULE common_dict

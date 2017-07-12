Attribute VB_Name = "callDLL"
Option Explicit

Public Declare PtrSafe Function LoadLibrary Lib "kernel32" Alias "LoadLibraryA" ( _
  ByVal lpLibFileName As String) As Long

Public Declare PtrSafe Function FreeLibrary Lib "kernel32" _
  (ByVal hLibModule As Long) As Long


Public Declare PtrSafe Function GetModuleHandle Lib "kernel32" Alias "GetModuleHandleA" ( _
  ByVal lpModuleName As String) As Long

'exp

Public Declare PtrSafe Function srmdllExpEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllExpMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllExpRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllExpInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllExpMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'gamma

Public Declare PtrSafe Function srmdllGammaEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, ByRef divide As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllGammaMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllGammaRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllGammaInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllGammaMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'pareto

Public Declare PtrSafe Function srmdllParetoEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllParetoMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllParetoRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllParetoInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllParetoMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'truncated normal

Public Declare PtrSafe Function srmdllTNormEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllTNormMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTNormRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTNormInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTNormMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'lognormal

Public Declare PtrSafe Function srmdllLNormEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllLNormMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLNormRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLNormInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLNormMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'truncated logistic

Public Declare PtrSafe Function srmdllTLogistEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllTLogistMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTLogistRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTLogistInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTLogistMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'log logistic

Public Declare PtrSafe Function srmdllLLogistEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllLLogistMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLLogistRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLLogistInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLLogistMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'truncated extreme value max

Public Declare PtrSafe Function srmdllTXvMaxEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllTXvMaxMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTXvMaxRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTXvMaxInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTXvMaxMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'log extreme value max

Public Declare PtrSafe Function srmdllLXvMaxEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllLXvMaxMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLXvMaxRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLXvMaxInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLXvMaxMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'truncated extreme value min

Public Declare PtrSafe Function srmdllTXvMinEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllTXvMinMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTXvMinRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTXvMinInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllTXvMinMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

'log extreme value min

Public Declare PtrSafe Function srmdllLXvMinEMstep _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef num As Double, ByRef ftype As Long, _
   ByRef npara As Long, ByRef para As Double, ByRef diff As Double, ByRef llf As Double, ByRef total As Double) As Long

Public Declare PtrSafe Function srmdllLXvMinMVF _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef mean As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLXvMinRate _
  Lib "srats2010.dll" _
  (ByRef dsize As Long, ByRef time As Double, ByRef rate As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLXvMinInverseMVF _
  Lib "srats2010.dll" _
  (ByRef value As Double, ByRef time As Double, _
   ByRef npara As Long, ByRef para As Double) As Long

Public Declare PtrSafe Function srmdllLXvMinMTTF _
  Lib "srats2010.dll" _
  (ByRef stime As Double, ByRef h As Double, ByRef m As Double, _
   ByRef npara As Long, ByRef para As Double) As Long






//Copyright 2019 Andrey S. Ionisyan (anserion@gmail.com)
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

//=====================================================================
//арифметика кватернионов
//=====================================================================

unit arith_quaternion;
{$mode objfpc}{$H+}
interface
uses math;
type 
TQuaternion=record
   a,b,c,d:real;
end;

TQuaternionVector=array of TQuaternion;
TQuaternionMatrix=array of TQuaternionVector;

function q_zero:TQuaternion;
procedure q_set(var q:TQuaternion; a,b,c,d:real);
function q_rnd:TQuaternion;
//function q_root_of_one():TQuaternion;
function q_dup(q:TQuaternion):TQuaternion;
function q_amp(q:TQuaternion):real;
function q_amp_bcd(q:TQuaternion):real;
function q_arg(q:TQuaternion):real;
function q_norm(q:TQuaternion):TQuaternion;

function q_conj(q:TQuaternion):TQuaternion;
function q_neg(q:TQuaternion):TQuaternion;
function q_inv(q:TQuaternion):TQuaternion;

function q_add(a,b:TQuaternion):TQuaternion;
function q_sub(a,b:TQuaternion):TQuaternion;

function q_scalar(lambda:real; q:TQuaternion):TQuaternion;
function q_dot(a,b:TQuaternion):real;
function q_dot_q(a,b:TQuaternion):TQuaternion;
function q_dot_bcd(a,b:TQuaternion):real;
function q_cross(a,b:TQuaternion):TQuaternion;
function q_cross_bcd(a,b:TQuaternion):TQuaternion;
function q_outer(a,b:TQuaternion):TQuaternion;
function q_mul(a,b:TQuaternion):TQuaternion;
function q_div(a,b:TQuaternion):TQuaternion;

function q_cos(q1,q2:TQuaternion):real;
function q_i(q:TQuaternion):TQuaternion;
function q_exp(q:TQuaternion):TQuaternion;
function q_ln(q:TQuaternion):TQuaternion;

implementation

procedure q_set(var q:TQuaternion; a,b,c,d:real);
begin q.a:=a; q.b:=b; q.c:=c; q.d:=d; end;

function q_zero:TQuaternion;
begin q_zero.a:=0; q_zero.b:=0; q_zero.c:=0; q_zero.d:=0; end;

function q_rnd:TQuaternion;
begin q_set(q_rnd,random,random,random,random); end;

function q_dup(q:TQuaternion):TQuaternion;
begin q_dup.a:=q.a; q_dup.b:=q.b; q_dup.c:=q.c; q_dup.d:=q.d; end;

function q_conj(q:TQuaternion):TQuaternion;
begin q_conj.a:=q.a; q_conj.b:=-q.b; q_conj.c:=-q.c; q_conj.d:=-q.d; end;

function q_neg(q:TQuaternion):TQuaternion;
begin q_neg.a:=-q.a; q_neg.b:=-q.b; q_neg.c:=-q.c; q_neg.d:=-q.d; end;

function q_inv(q:TQuaternion):TQuaternion;
begin
   if (q.a=0) and (q.b=0) and (q.c=0) and (q.d=0)
   then q_inv:=q_zero
   else q_inv:=q_scalar(1/(q.a*q.a+q.b*q.b+q.c*q.c+q.d*q.d),q_conj(q));
end;

function q_amp(q:TQuaternion):real;
begin q_amp:=sqrt(q.a*q.a+q.b*q.b+q.c*q.c+q.d*q.d); end;

function q_amp_bcd(q:TQuaternion):real;
begin q_amp_bcd:=sqrt(q.b*q.b+q.c*q.c+q.d*q.d); end;

function q_arg(q:TQuaternion):real;
begin
   if (q.a=0) and (q.b=0) and (q.c=0) and (q.d=0)
   then q_arg:=0
   else q_arg:=arccos(q.a/q_amp(q));
end;

function q_norm(q:TQuaternion):TQuaternion;
begin
   if (q.a=0) and (q.b=0) and (q.c=0) and (q.d=0)
   then q_norm:=q
   else q_norm:=q_scalar(1/q_amp(q),q);
end;

function q_add(a,b:TQuaternion):TQuaternion;
begin
   q_add.a:=a.a+b.a;
   q_add.b:=a.b+b.b;
   q_add.c:=a.c+b.c;
   q_add.d:=a.d+b.d;
end;

function q_sub(a,b:TQuaternion):TQuaternion;
begin
   q_sub.a:=a.a-b.a;
   q_sub.b:=a.b-b.b;
   q_sub.c:=a.c-b.c;
   q_sub.d:=a.d-b.d;
end;

function q_scalar(lambda:real; q:TQuaternion):TQuaternion;
begin
   q_scalar.a:=lambda*q.a;
   q_scalar.b:=lambda*q.b;
   q_scalar.c:=lambda*q.c;
   q_scalar.d:=lambda*q.d;
end;

function q_dot(a,b:TQuaternion):real;
begin q_dot:=a.a*b.a+a.b*b.b+a.c*b.c+a.d*b.d; end;

function q_dot_bcd(a,b:TQuaternion):real;
begin q_dot_bcd:=a.b*b.b+a.c*b.c+a.d*b.d; end;

function q_cross(a,b:TQuaternion):TQuaternion;
begin q_cross:=q_scalar(0.5,q_sub(q_mul(a,b),q_mul(b,a))); end;

function q_dot_q(a,b:TQuaternion):TQuaternion;
begin
   q_dot_q:=q_scalar(0.5,q_add(q_mul(q_conj(a),b),q_mul(q_conj(b),a)));
end;

function q_cross_bcd(a,b:TQuaternion):TQuaternion;
begin
   q_cross_bcd.a:=0;
   q_cross_bcd.b:=a.c*b.d-a.d*b.c;
   q_cross_bcd.c:=a.d*b.b-a.b*b.d;
   q_cross_bcd.d:=a.b*b.c-a.c*b.b;
end;

function q_outer(a,b:TQuaternion):TQuaternion;
begin
   q_outer:=q_scalar(0.5,q_sub(q_mul(q_conj(a),b),q_mul(q_conj(b),a)));
end;

function q_mul(a,b:TQuaternion):TQuaternion;
var tmp1,tmp2,tmp3:TQuaternion;
begin
   q_mul.a:=a.a*b.a-q_dot_bcd(a,b);
   
   tmp1:=q_scalar(a.a,b);
   tmp2:=q_scalar(b.a,a);
   tmp3:=q_cross_bcd(a,b);
   
   q_mul.b:=tmp1.b+tmp2.b+tmp3.b;
   q_mul.c:=tmp1.c+tmp2.c+tmp3.c;
   q_mul.d:=tmp1.d+tmp2.d+tmp3.d;
end;

function q_div(a,b:TQuaternion):TQuaternion;
begin q_div:=q_mul(a,q_inv(b)); end;

function q_cos(q1,q2:TQuaternion):real;
var amp1,amp2:real;
begin
   amp1:=q_amp(q1); amp2:=q_amp(q2);
   if (amp1=0) or (amp2=0)
   then q_cos:=0
   else q_cos:=q_dot(q1,q2)/(amp1*amp2);
end;

function q_i(q:TQuaternion):TQuaternion;
begin
   if (q.b=0) and (q.c=0) and (q.d=0)
   then q_i:=q
   else
      begin
         q.a:=0;
         q_i:=q_scalar(1.0/q_amp(q),q);
      end;
end;

function q_exp(q:TQuaternion):TQuaternion;
var amp_u:real; tmp,u,i:TQuaternion;
begin
   if (q.a=0) and (q.b=0) and (q.c=0) and (q.d=0)
   then begin tmp.a:=1; tmp.b:=0; tmp.c:=0; tmp.d:=0; q_exp:=tmp; end
   else 
      begin
         u:=q_dup(q); u.a:=0;
         amp_u:=q_amp(u);
         i:=q_i(q);
         tmp:=q_scalar(sin(amp_u),i);
         tmp.a:=tmp.a+cos(amp_u);
         q_exp:=q_scalar(exp(amp_u),tmp);
      end;
end;

function q_ln(q:TQuaternion):TQuaternion;
var amp_q:real; tmp:TQuaternion;
begin
   if (q.a=0) and (q.b=0) and (q.c=0) and (q.d=0)
   then q_ln:=q_zero
   else
      begin
         amp_q:=q_amp(q);
         tmp:=q_scalar(q_arg(q),q_i(q));
         tmp.a:=tmp.a+ln(amp_q);
         q_ln:=tmp;
      end;
end;

end.

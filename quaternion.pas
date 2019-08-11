program quaternion;
uses arith_quaternion;

procedure q_print(s_pre:string; q:TQuaternion; s_post:string);
begin
   write(s_pre);
   write(q.a:8:3,q.b:8:3,q.c:8:3,q.d:8:3);
   writeln(s_post);
end;

var a,b,c,u,i,i2:TQuaternion;
begin
   q_set(a,2,1,0,0); q_set(b,3,0,1,0);
   
   a:=q_rnd; b:=q_rnd; a.a:=0; b.a:=0;
   q_print('a=(',a,')');
   q_print('b=(',b,')');
   writeln('amp_a=',q_amp(a):8:3);
   writeln('amp_b=',q_amp(b):8:3);

   c:=q_norm(a);
   q_print('norm(a)=(',c,')');
   writeln('amp_c=',q_amp(c):8:3);

   c:=q_norm(b);
   q_print('norm(b)=(',c,')');
   writeln('amp_c=',q_amp(c):8:3);
   
   c:=q_mul(a,b);
   q_print('a*b=(',c,')');
   writeln('amp_c=',q_amp(c):8:3);

   c:=q_div(a,b);
   q_print('a/b=(',c,')');
   writeln('amp_c=',q_amp(c):8:3);

   c:=q_dot_q(a,b);
   q_print('a dot b=(',c,')');

   c:=q_cross(a,b);
   q_print('a cross b=(',c,')');

   u:=q_dup(a); u.a:=0;
   i:=q_scalar(1/q_amp(u),u);
   q_print('!u!^(-1) * u = (',i,')');
   
   i2:=q_mul(i,i);
   q_print('i^2 = (',i2,')');
   
   q_set(a,1,0,0,0); q_print('a=(',a,')');
   q_print('exp(a) = (',q_exp(a),')');
   q_print('ln(a) = (',q_ln(a),')');
   
   q_print('exp(ln(a)) = (',q_exp(q_ln(a)),')');
end.

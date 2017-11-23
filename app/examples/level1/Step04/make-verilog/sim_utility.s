
global super_start :label;
global stop :label;

extern start :label;

begin ".load"

    end ".load";
        
begin ".init"
<super_start>
    goto start;
<stop>
    [1003_3fffh] = gr7;
    halt;
global _init: label;
<_init>
end ".init";

begin ".fini"
global _fini: label;
<_fini>
end ".fini";



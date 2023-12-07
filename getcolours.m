function [colour]=getcolours(ind)
    ind = 1+ mod(ind-1,24);

    ClassicBlue = [52,86,139];
    LivingCoral = [255,111,97];
    Greenery = [136,176,75];
    RoseQuartz = [247,202,201];
    Serenity = [146,168,209];

    Marsala = [149,82,81];
    UltraViolet = [107,91,149];
    RadiandOrchid = [181,101,167];
    Emerald = [0,155,119];
    TangerineTango = [221,65,36];

    Honeysucle = [214,80,118];
    Turquiose = [68,184,172];
    Mimosa = [239,192,80];
    BlueIzis = [91,94,166];
    ChiliPepper = [155,35,53];

    SandDollar = [223,207,190];
    BlueTurquoise = [85,180,176];
    Tigerlily = [225,93,68];
    AguaSky = [127,205,205];
    TrueRed = [188,36,60];

    FuchsiaRose = [152,180,212];
    Marigold = [253,172,83];
    Cerulean = [155,183,212];
    Rust = [181,90,48];
    FrenchBlue = [0,114,181];


    colours = [ClassicBlue; ...
      LivingCoral;  Greenery;...
      RoseQuartz;  UltraViolet;...
      SandDollar; Serenity ;        ...
      RadiandOrchid ;   Marigold;  Emerald ;     ...
      Honeysucle ;...
      Turquiose ;        BlueIzis ; Cerulean; Tigerlily;   ChiliPepper ;     ...
      BlueTurquoise ;  Mimosa ;     AguaSky;  TangerineTango ; Marsala;    TrueRed;    ...
        Rust; FrenchBlue];

    colour = colours(ind,:)/256;
end
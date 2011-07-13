% initialized variables

function [Gvars input] = setupVars()

Gvars = config;
input = userInput(Gvars);
[Gvars, input] = initControlVars(Gvars, input);

end
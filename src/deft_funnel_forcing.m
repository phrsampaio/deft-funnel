function value = deft_funnel_forcing(number, argument)

  if number == 1
    value = 1.0e-2*min(1.0, min(abs(argument), argument^2));
  elseif number == 2
    value = 1.0e-2*min(1.0, min(abs(argument), argument^2)); % TDB
  elseif number == 3
    value = 1.0e-2*min(1.0, min(abs(argument), argument^2)); % TDB
  elseif number == 4
    value = min(1.0, argument);
  end
  
end % end of deft_funnel_forcing

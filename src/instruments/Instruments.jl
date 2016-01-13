## TOP LEVEL CALCULATION METHODS - KEEPING TRACK OF CALCULATION STATE ##
function calculate!{I <: Instrument}(inst::I)
  if !inst.calculated
    inst.calculated = true
    perform_calculations!(inst)
  end

  return inst
end

function recalculate!{I <: Instrument}(inst::I)
  inst.calculated = false
  calculate!(inst)

  return inst
end

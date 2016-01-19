## TOP LEVEL CALCULATION METHODS - KEEPING TRACK OF CALCULATION STATE ##
function update_pricing_engine!{I <: Instrument, P <: PricingEngine}(inst::I, pe::P)
  inst.pricingEngine = pe
  inst.lazyMixin.calculated = false
  return inst
end

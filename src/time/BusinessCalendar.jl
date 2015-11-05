# Business Calendars (adapted from Ito.jl and BusinessDays.jl)

abstract BusinessCalendar

abstract WesternCalendar <: BusinessCalendar
abstract OrthodoxCalendar <: BusinessCalendar

# US Calendars
abstract UnitedStatesCalendar <: WesternCalendar

type USSettlementCalendar <: UnitedStatesCalendar; end
type USNYSECalendar <: UnitedStatesCalendar; end
type USNERCCalendar <: UnitedStatesCalendar; end
type USGovernmentBondCalendar <: UnitedStatesCalendar; end

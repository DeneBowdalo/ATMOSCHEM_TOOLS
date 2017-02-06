Program shuffle
!
! Shuffle a deck of 52 cards and generate 4 hands.
! The method is to pick a card at random and store it.
! Then that card is moved to the high end of the list
! So when 2 cards have been used cards 1-50 will be the
! ones we can choose from, cards 51, 52 will be the
! ones we have already used.
!
      IMPLICIT NONE
      INTEGER ::i, card_selected, selected_value, selected_suit
      INTEGER :: hand
      INTEGER, DIMENSION(52) :: cards
      REAL :: a_number
      CHARACTER(len=1),DIMENSION(4)  :: suit
      CHARACTER(len=1),DIMENSION(13) :: value

!
! Use implied do loops to assign vaules to the arrays
!
      cards = (/ (i,i=1,52) /)
      suit = (/ 'C','D','H','S' /)
      value = (/ '1','2','3','4','5','6','7','8','9','T','J','Q','K' /)
      hand = 0
!
! Loop through the 52 cards
! Every 14th card move onto a new hand
!
      DO i = 52,2,-1
         IF(MOD(i,13).eq.0)THEN
            hand = hand + 1
            WRITE(*,*)'Hand',hand
            WRITE(*,*)'***************'
         END IF
!
! Generate a random number between 0 and 1
! and use it to pick a suit and value
!
         CALL random_number(a_number)
         card_selected = INT(a_number*REAL(i))+1
         selected_suit = (cards(card_selected)-1)/13 + 1
         selected_value = MOD(cards(card_selected)-1,13) + 1
         WRITE(*,*)value(selected_value),suit(selected_suit)
!
! Move the card into the list of cards already used
!
         IF(card_selected.ne.i)THEN
             cards(card_selected) = cards(i)
         END IF
      END DO
!
! Assign the remaining card
!
      selected_suit = (cards(1)-1)/13 + 1
      selected_value = MOD(cards(1)-1,13) + 1
      WRITE(*,*)value(selected_value),suit(selected_suit)
 END PROGRAM shuffle
